/* (C) 2013-2014 TU Delft
   (C) 2013-2014 AHU Holding BV */

#define __STDC_FORMAT_MACROS
#include <tclap/CmdLine.h>
#include <stdio.h>
#include <iostream>
#include <stdint.h>
#include <stdlib.h>
#include <fstream>
#include <map>
#include <unordered_map>
#include <boost/math/distributions/poisson.hpp>
#include <vector>
#include <string>
#include <string.h>
#include <stdexcept>
#include <forward_list>
#include <inttypes.h>
#include <algorithm>
#include <numeric>

#include <errno.h>
#include <math.h>
#include <signal.h>
#include <boost/progress.hpp>
#include <boost/format.hpp>
#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include "dnamisc.hh"
#include <fenv.h>
#include <memory>
#include <sstream>
#include "zstuff.hh"
#include "geneannotated.hh"
#include "misc.hh"
#include "fastq.hh"
#include <mba/diff.h>
#include <mba/msgno.h>
#include "antonie.hh"
#include "saminfra.hh"
#include "refgenome.hh"
#include "compat.hh"

extern "C" {
#include "hash.h"
}
using boost::lexical_cast;

/*! \mainpage Antonie DNA Software
  The Antonie DNA software reads DNA reads as FastQRead objects, and matches them to a
  ReferenceChromosome. The reading happens through a StereoFASTQReader instance, which in turn
  has two FASTQReader instances, thus supporting paired end reads.

  Successfully matched reads are written out through the SAMWriter class, if so configured.
*/

using namespace std;
using namespace boost::algorithm;
namespace io = boost::iostreams;
typedef io::tee_device<std::ostream, std::ostringstream> TeeDevice;
typedef io::stream< TeeDevice > TeeStream;
TeeStream* g_log;
string g_name;

void ReferenceChromosome::printCoverage(FILE* jsfp, const std::string& histoName)
{
  uint64_t totCoverage=0, noCoverages=0;
  unsigned int cov;

  bool wasNul=true;
  string::size_type prevNulpos=0;

  vector<unsigned int> covhisto;
  covhisto.resize(1000);
  VarMeanEstimator vmeDepth;
  d_unmRegions.clear();
  //  ofstream covfile("coverage");
  //  int prevcov = d_mapping[0].coverage;

  map<double, VarMeanEstimator> gcCoverage;

  for(string::size_type pos = 0; pos < d_mapping.size(); ++pos) {
    cov = d_mapping[pos].coverage;
    string gcsnip=snippet((pos > 20) ? (pos - 20) : 1, pos+20);
    double gc=getGCContent(gcsnip); // 1 based!
    
    if(gcsnip.length()==40) {// we get strange results otherwise
      gcCoverage[gc](cov);
    }
    else {
      // cout <<"Odd gcsnip length: "<<gcsnip.length()<<", pos = "<<pos<<endl;
      
    }

    //    covfile << pos << '\t' << cov << '\t' << ((int)cov)-prevcov << '\t' << gc  <<'\t' << getGCContent(snippet(pos > 7 ? pos-7 : 1, pos+8)) << '\n';
    // prevcov=cov;

    vmeDepth(cov);

    bool noCov = cov < 5;
    if(cov >= covhisto.size()) {
      covhisto.resize(2*cov+1);
    }

    covhisto[cov]++;
    
    totCoverage += cov;

    if(noCov) {
      noCoverages++;
    }
    
    if(!noCov && wasNul) {
      if(prevNulpos > 40 && pos + 40 < d_genome.length()) {
	Unmatched unm;
	unm.left = d_genome.substr(prevNulpos-40, 40);
	unm.right = d_genome.substr(pos, 40);
	unm.unmatched = d_genome.substr(prevNulpos, pos-prevNulpos);
	unm.pos = prevNulpos;
	d_unmRegions.push_back(unm);
      }
      wasNul=false;
    }
    else if(noCov && !wasNul) {
      wasNul=true;
      prevNulpos = pos;
    }
  }

  ofstream gcv("gccoverage");
  vector<pair<double, double> > gc, gclo, gchi;
  for(const auto& ent : gcCoverage) {
    gc.push_back({ent.first, mean(ent.second)});
    gclo.push_back({ent.first, mean(ent.second) - sqrt(variance(ent.second))});
    gchi.push_back({ent.first, mean(ent.second) + sqrt(variance(ent.second))});
    gcv << ent.first << '\t' << mean(ent.second) << '\t' << sqrt(variance(ent.second))<<endl;
  }
  fprintf(jsfp, "var gccov=%s;\nvar gccovlo=%s;\nvar gccovhi=%s;\n", jsonVectorPair(gc).c_str(), jsonVectorPair(gclo).c_str(), jsonVectorPair(gchi).c_str());
  
  Clusterer<Unmatched> cl(100);
  
  for(auto unm : d_unmRegions) {
    cl.feed(unm);
  }

  (*g_log) << (boost::format("Average depth: %|40t|    %10.2f +- %.2f\n") % mean(vmeDepth) % sqrt(variance(vmeDepth))).str();

  double expected=-1; // size()*(1-erf(mean(vmeDepth)/(sqrt(variance(vmeDepth))*sqrt(2))));
  (*g_log) << (boost::format("Undercovered nucleotides: %|40t| %10d (%.2f%%), %d ranges, %.1f could be expected\n") % noCoverages % (noCoverages*100.0/d_mapping.size()) % cl.d_clusters.size() %expected).str();
  
  uint64_t total = std::accumulate(covhisto.begin(), covhisto.end(), 0), cumul=0;

  // snip off once we have 99.9%
  for(auto iter = covhisto.begin(); iter != covhisto.end(); ++iter) {
    cumul += *iter;
    if(cumul > total*0.999) {
      covhisto.resize((iter - covhisto.begin()));
      break;
    }
  }

  fputs(jsonVector(covhisto, histoName, [&total](dnapos_t dp) { return 1.0*dp/total; }).c_str(), jsfp);
}


//! 0 if nothing interesting, positive if our read has insert at that position, negative if we have a delete at that position
int MBADiff(dnapos_t pos, const FastQRead& fqr, const string& reference, unsigned int* amount=0)
{
  string::size_type n, m;
  int d,  sn;
  struct varray *ses = varray_new(sizeof(struct diff_edit), NULL);
  
  n = reference.length();
  m = fqr.d_nucleotides.length();
  if ((d = diff(fqr.d_nucleotides.c_str(), 0, m, reference.c_str(), 0, n, NULL, NULL, NULL, 0, ses, &sn, NULL)) == -1) {
    MMNO(errno);
    printf("Error\n");
    return EXIT_FAILURE;
  }
  int ret=0;

  if(sn > 3 && sn < 10) {
    auto *match1 = (struct diff_edit*)varray_get(ses, 0),
      *change1=(struct diff_edit*)varray_get(ses, 1), 
      *match2=(struct diff_edit*)varray_get(ses, 2);

    int totMatch = match1->len + match2->len;
    if(match1->op == DIFF_MATCH && totMatch > 10*change1->len &&  match1->len > 10 && match2->len > 10 &&
       match2->op == DIFF_MATCH && (d - 2*change1->len) < 4) {
      if(0 && amount )
		printf("pos %u, d=%u (%d) sn=%d, rpos=%u, calling '%s' of %d bases\nUS:  %s\nREF: %s\n", pos, d, d - 2*change1->len, 
	     sn, pos+change1->off, change1->op==DIFF_INSERT ? "Delete" : "Insert", change1->len, fqr.d_nucleotides.c_str(), reference.c_str());
      
      if(amount)
	*amount=change1->len;

      if(change1->op == DIFF_INSERT) 
	ret=-change1->off;
      else
	ret=change1->off;
#if 0
      for (int i = 0; i < sn; i++) {
	struct diff_edit *e = (struct diff_edit*)varray_get(ses, i);
	
	switch (e->op) {
	case DIFF_MATCH:
	  printf("MAT: ");
	  fwrite(fqr.d_nucleotides.c_str() + e->off, 1, e->len, stdout);
	  break;
	case DIFF_INSERT:
	  printf("INS: ");
	  fwrite(reference.c_str() + e->off, 1, e->len, stdout);
	  break;
	case DIFF_DELETE:
	  printf("DEL: ");
	  fwrite(fqr.d_nucleotides.c_str() + e->off, 1, e->len, stdout);
	  break;
	}
	printf("\n");
      }
#endif
    }
  }

  varray_del(ses);

  return ret;
}

unsigned int diffScore(ReferenceChromosome& rg, dnapos_t pos, const FastQRead& fqfrag, int qlimit)
{
  unsigned int diffcount=0;
  string reference = rg.snippet(pos, pos + fqfrag.d_nucleotides.length());
  for(string::size_type i = 0; i < fqfrag.d_nucleotides.size() && i < reference.size();++i) {
    if(fqfrag.d_nucleotides[i] != reference[i] && fqfrag.d_quality[i] > qlimit) 
      diffcount++;
  }

  if(diffcount >= 5) { // bit too different, try mbadiff!
    int res=MBADiff(pos, fqfrag, reference);
    if(res < 0 || res > 0)
      return 1;
  }

  return diffcount;
}

vector<dnapos_t> getTriplets(const vector<pair<dnapos_t, char>>& together, unsigned int interval, unsigned int shift) 
{
  vector<dnapos_t> ret;

  for(unsigned int i = 0; i < together.size() - 2; ++i) {
    if(together[i].second=='L' && together[i+1].second=='M' && together[i+2].second=='R' && 
       together[i+1].first - together[i].first < 1.2*interval && together[i+2].first - together[i+1].first < 1.2*interval) {
      dnapos_t lpos;
      lpos=together[i].first; 

      if(lpos < shift)
        continue;
      
      ret.push_back(lpos-shift);
    }
  }         
  return ret;
}

void printCorrectMappings(FILE* jsfp, const ReferenceChromosome& rg, const std::string& name)
{
  fprintf(jsfp, "%s=[", name.c_str());
  for(unsigned int i=0; i < rg.d_correctMappings.size() ;++i) {
    if(!rg.d_correctMappings[i] || !rg.d_wrongMappings[i])
      continue;
    double total=rg.d_correctMappings[i] + rg.d_wrongMappings[i];
    double error= rg.d_wrongMappings[i]/total;
    double qscore=-10*log10(error);
    fprintf(jsfp, "%s[%d,%.2f]", i ? "," : "", i, qscore);
    //    cout<<"total "<<total<<", error: "<<error<<", qscore: "<<qscore<<endl;
  }
  fprintf(jsfp,"];\n");
}


//! Keeps a tally of correct and incorrect mappings
struct qtally
{
  qtally() : correct{0}, incorrect{0}{}
  uint64_t correct;
  uint64_t incorrect;
};

int MapToReference(ReferenceChromosome& rg, dnapos_t pos, FastQRead fqfrag, int qlimit, BAMWriter* sbw, vector<qtally>* qqcounts, int* outIndel=0)
{
  if(pos > rg.size()) // can happen because of inserts or circular genomes
    return false;
  if(outIndel)
    *outIndel=0;
  string reference = rg.snippet(pos, pos + fqfrag.d_nucleotides.length());

  double diffcount=0;
  for(string::size_type i = 0; i < fqfrag.d_nucleotides.size() && i < reference.size();++i) {
    if(fqfrag.d_nucleotides[i] != reference[i]) {
      if(fqfrag.d_quality[i] > qlimit) 
	diffcount++;
      else
	diffcount+=0.5;
    }
  }
  bool didMap=false;

  if(diffcount < 5) {
    didMap=true;
    rg.mapFastQ(pos, fqfrag);
    if(sbw)
      sbw->qwrite(pos, fqfrag);
  }
  else {
    unsigned int amount=0;
    int indel=MBADiff(pos, fqfrag, reference, &amount);
    if(outIndel)
      *outIndel=indel;

    if(indel) {
      rg.mapFastQ(pos, fqfrag, indel);
      if(sbw)
	sbw->qwrite(pos, fqfrag, indel);
      didMap=true;

      diffcount=1; // makes sure we get mapped anyhow
      if(indel > 0) { // our read has an insert at this position
	//	cout<<"Mapping an insert: "<<pos+indel<<" of "<<amount<<" codons, "<<fqfrag.d_nucleotides.substr(indel, amount)<<endl;

	// down below, everything will match, so we need to locimap here
        rg.d_locimap[pos+indel].samples.push_back({fqfrag.d_nucleotides[indel], fqfrag.d_quality[indel], 
	      (bool)(fqfrag.reversed ^ ((unsigned int)indel > fqfrag.d_nucleotides.length()/2)),      // head or tail
	      fqfrag.d_nucleotides.substr(indel, amount)});

        fqfrag.d_nucleotides.erase(indel, amount); // this makes things align again
        fqfrag.d_quality.erase(indel, amount); 
        rg.d_insertCounts[pos+indel]++; 
      } else {      // our read has an erase at this position
	//	cout<<"Mapping a delete at "<<pos-indel<<" of " <<amount<<" codons"<<endl;
        fqfrag.d_nucleotides.insert(-indel, amount, 'X');
	//	cout<<"REF: "<<reference<<endl<<"US:  "<<fqfrag.d_nucleotides<<endl;
        fqfrag.d_quality.insert(-indel, amount, 40);
      }
    }
  }

  unsigned int readMapPos;
  for(string::size_type i = 0; i < fqfrag.d_nucleotides.size() && i < reference.size();++i) {
    readMapPos = fqfrag.reversed ? ((reference.length()- 1) - i) : i; // d_nucleotides might have an insert
      
    char c =  fqfrag.d_nucleotides[i];

    if(c != reference[i]) {
      //      diff.append(1, fqfrag.d_quality[i] > qlimit ? '!' : '^');
      //      cout<<"Have diff at "<<pos+i<<", diffCount="<<diffcount<<", qfilt="<< (fqfrag.d_quality[i] > qlimit)<<endl;
      if(fqfrag.d_quality[i] > qlimit && diffcount < 5) 
        rg.d_locimap[pos+i].samples.push_back({fqfrag.d_nucleotides[i], fqfrag.d_quality[i], 
	      (bool)(fqfrag.reversed ^ (i > fqfrag.d_nucleotides.length()/2))}); // head or tail
      
      if(diffcount < 5) {
	unsigned int q = (unsigned int)fqfrag.d_quality[i];
	(*qqcounts)[q].incorrect++;
	rg.d_wrongMappings[readMapPos]++;
      }
    }
    else {
      // diff.append(1, ' ');
      rg.cover(pos+i,fqfrag.d_quality[i], qlimit);
      if(diffcount < 5) {
	(*qqcounts)[(unsigned int)fqfrag.d_quality[i]].correct++;
	rg.d_correctMappings[readMapPos]++;
      }
    }
  }
  //  if(diffcount > 5) {
  //  cout<<"US:  "<<fqfrag.d_nucleotides<<endl<<"DIF: ";
  //  cout<<diff<<endl<<"REF: "<<reference<<endl;
  //  cout<<"QUA: "<<fqfrag.d_quality<<endl;
  //  cout<<endl;
  //}
  //return diff;
  return didMap;
}


void emitRegion(FILE*fp, ReferenceChromosome& rg, StereoFASTQReader& fastq, const string& name, unsigned int index, dnapos_t start, 
		dnapos_t stop, const std::string& report_="", int maxVarcount=-1)
{
  if(stop > rg.size())
    stop = rg.size();
  if(start > rg.size())
    start=1;
  dnapos_t dnapos = (start+stop)/2;
  fprintf(fp, "region[%d]={reference: '%s', name:'%s', pos: %d, depth: [",  index, rg.d_name.c_str(), name.c_str(), dnapos);
  for(dnapos_t pos = start; pos < stop; ++pos) {
    if(pos != start) 
      fprintf(fp, ",");
    fprintf(fp, "[%d,%d]", pos, rg.d_mapping[pos].coverage);
  }
  fprintf(fp, "], ");
  vector<double> aProb(stop-start), cProb(stop-start), gProb(stop-start), tProb(stop-start), xProb(stop-start);
  for(dnapos_t pos = start; pos < stop; ++pos) {
    if(rg.d_locimap.count(pos)) {
      for(const auto& locus: rg.d_locimap[pos].samples) {
	acgtxDo(locus.nucleotide, 
		[&]() { aProb[pos-start]+= locus.quality/30.0; }, 
		[&]() { cProb[pos-start]+= locus.quality/30.0; }, 
		[&]() { gProb[pos-start]+= locus.quality/30.0; }, 
		[&]() { tProb[pos-start]+= locus.quality/30.0; },
		[&]() { xProb[pos-start]+= locus.quality/30.0; });
      }
    }
  }
  
  fprintf(fp, "aProb: %s, cProb: %s, gProb: %s, tProb: %s, xProb: %s",
	  jsonVectorX(aProb, [start](int i){return i+start;}).c_str(), 
	  jsonVectorX(cProb, [start](int i){return i+start;}).c_str(), 
	  jsonVectorX(gProb, [start](int i){return i+start;}).c_str(), 
	  jsonVectorX(tProb, [start](int i){return i+start;}).c_str(),
	  jsonVectorX(xProb, [start](int i){return i+start;}).c_str());

  string picture; // =rg.getMatchingFastQs(start, stop, fastq);
  string snippet=rg.snippet(start, dnapos) + " | " +rg.snippet(dnapos, stop);
  replace_all(picture, "\n", "\\n");
  string report = replace_all_copy(report_, "\n", "\\n");

  string annotations;
  int gene=0;
  if(rg.d_gar) {
    auto gas=rg.d_gar->lookup("", dnapos);
    abort(); // instead of "", it needs to have a name of a chromosome
    for(auto ga : gas) {
      replace_all(ga.tag, "'", "\\'");
      annotations += ga.name+" [" + ga.tag  + "], ";
      if(ga.gene)
	gene=1;
    }
  }
  replace_all(report, "'", "\\'");
  fprintf(fp,"picture: '%s', snippet: '%s', maxVarcount: %d, gene: %d, annotations: '%s', report: '%s'};\n", "", snippet.c_str(), maxVarcount, gene, annotations.c_str(), report.c_str());
  
  fputs("\n", fp);
  fflush(fp);
}

void emitRegion(FILE*fp, ReferenceChromosome& rg, StereoFASTQReader& fastq, const string& name, unsigned int index, dnapos_t start, const std::string& report="")
{
  emitRegion(fp, rg, fastq, name, index, start > 200 ? start-200 : 1, (start +200) < rg.size() ? (start + 200) : rg.size(), report);
}

unsigned int variabilityCount(const ReferenceChromosome& rg, dnapos_t position, const ReferenceChromosome::LociStats& lc, double* fraction)
{
  vector<int> counts(256);
  counts[rg.snippet(position, position+1)[0]]+=rg.d_mapping[position].coverage;
  
  int forwardCount=0;

  for(const auto& j : lc.samples) {
    counts[j.nucleotide]++;
    if(j.headOrTail)
      forwardCount++;
  }
  sort(counts.begin(), counts.end());
  unsigned int nonDom=0;
  for(unsigned int i=0; i < 255; ++i) {
    nonDom+=counts[i];
  }
  
  if(nonDom + counts[255] < 20) // depth
    return 0;

  *fraction = 1.0*forwardCount / (1.0*lc.samples.size());
  if(*fraction < 0.05 || *fraction > 0.95)
    return 0;

  return 100*nonDom/(nonDom+counts[255]);
}


vector<ReferenceChromosome::MatchDescriptor> fuzzyFind(FastQRead* fqfrag, ReferenceChromosome& rg, unsigned int keylen, int qlimit)
{
  vector<ReferenceChromosome::MatchDescriptor> ret;

  string left, middle, right;
  typedef pair<dnapos_t, char> tpos;
  vector<dnapos_t> lpositions, mpositions, rpositions;

  if(fqfrag->d_nucleotides.length() < 3*keylen) // too short
    return ret;

  unsigned int interval=(fqfrag->d_nucleotides.length() - 3*keylen)/3;

  for(unsigned int attempts=0; attempts < interval; attempts += 3) {
    for(int tries = 0; tries < 2; ++tries) {
      if(tries)
	fqfrag->reverse();
      left=fqfrag->d_nucleotides.substr(attempts, keylen);   
      middle=fqfrag->d_nucleotides.substr(interval+attempts, keylen); 
      right=fqfrag->d_nucleotides.substr(2*interval+attempts, keylen);
      lpositions=rg.getReadPositions(left); 

      if(lpositions.empty())
	continue;
      mpositions=rg.getReadPositions(middle); 
      if(mpositions.empty())
	continue;
      rpositions=rg.getReadPositions(right);
      
      if(lpositions.size() + mpositions.size() + rpositions.size() < 3)
	continue;
      vector<tpos> together;
      for(auto fpos: lpositions) { together.push_back(make_pair(fpos, 'L')); }        
      for(auto fpos: mpositions) { together.push_back(make_pair(fpos, 'M')); }
      for(auto fpos: rpositions) { together.push_back(make_pair(fpos, 'R')); }
      
      sort(together.begin(), together.end());
      
      auto matches=getTriplets(together, interval, attempts);
      //	random_shuffle(matches.begin(), matches.end()); // should prevent pileup because if 'score==0' shortcut below
      int score;
      for(auto match : matches) {
	if(std::find_if(ret.begin(), ret.end(), 
			[&match](const ReferenceChromosome::MatchDescriptor& md){ return md.pos==match;}) != ret.end())
	  continue;
	
	score = diffScore(rg, match, *fqfrag, qlimit);
        
	ret.push_back({&rg, match, fqfrag->reversed, score});
	if(score==0) // won't get any better than this
	  return ret;
      }
    }
  }
  return ret;
}

vector<ReferenceChromosome::MatchDescriptor> fuzzyFind(FastQRead* fqfrag, vector<unique_ptr<ReferenceChromosome> >& refs, int keylen, int qlimit)
{
  vector<ReferenceChromosome::MatchDescriptor> ret;
  for(auto& rg : refs) {
    auto inter = fuzzyFind(fqfrag, *rg, keylen, qlimit);
    for(auto& i : inter)
      ret.push_back(i); // XXX must be a better way
  }
  return ret;
}


typedef vector<VarMeanEstimator> qstats_t;

void writeUnmatchedReads(const vector<uint64_t>& unfoundReads, StereoFASTQReader& fastq)
{
  FILE *fp=fopen("unfound.fastq", "w");
  FastQRead fqfrag;
  for(const auto& pos :  unfoundReads) {
    fastq.getRead(pos, &fqfrag);
    fprintf(fp, "@%s\n%s\n+\n%s\n", fqfrag.d_header.c_str(), fqfrag.d_nucleotides.c_str(), fqfrag.getSangerQualityString().c_str());
  }
  fclose(fp);
}

struct ClusterLocus
{
  unsigned int pos;
  ReferenceChromosome::LociStats locistat;
  bool operator<(const ClusterLocus& rhs) const
  {
    return pos < rhs.pos;
  }
};

string makeAminoReport(ReferenceChromosome& rg, dnapos_t pos, const vector<GeneAnnotation>& gas, 
		       const ReferenceChromosome::LociStats& locistat, string* headline, string* body)
{ 
  string origCodon{"XXX"}, newCodon;
  string gene;
  int nucOffset=0;
  bool orfSense=0;
  int aminoNum=0;
  string fmt2("                  ");

  int aCount{0}, cCount{0}, gCount{0}, tCount{0};
  char c=rg.snippet(pos, pos+1)[0];
  acgtDo(c, 
	 [&](){aCount += rg.d_mapping[pos].coverage;},
	 [&](){cCount += rg.d_mapping[pos].coverage;},
	 [&](){gCount += rg.d_mapping[pos].coverage;},
	 [&](){tCount += rg.d_mapping[pos].coverage;}
	 );
  
  for(auto j = locistat.samples.begin(); 
      j != locistat.samples.end(); ++j) {
    c=j->nucleotide;
    acgtDo(c, [&](){ aCount++; }, [&](){ cCount++; }, [&](){ gCount++; }, [&](){ tCount++; } );
  }
  for(const auto& ga : gas) {
    if(ga.gene) {
      gene = rg.snippet(ga.startPos, ga.stopPos+1);
      orfSense = ga.strand;
      if(ga.strand) {
	aminoNum = (pos - ga.startPos) / 3;
	nucOffset = (pos - ga.startPos) % 3;
      }
      else {
	aminoNum = (ga.stopPos - pos) / 3;
	nucOffset = (ga.stopPos - pos) % 3;
	reverseNucleotides(&gene);
      }
      if((unsigned int)(aminoNum*3 + 3) < gene.size())
	origCodon = gene.substr(aminoNum*3, 3);
    }
  }
  ostringstream ret;
  string residueString = lexical_cast<string>(1+aminoNum) + '/'+lexical_cast<string>(gene.size()/3);
  ret<<"Original codon: "<<origCodon<<", amino acid: "<<AminoAcidName(DNAToAminoAcid(origCodon.c_str()))<<", Residue "<<residueString<<", offset in codon "<<nucOffset<<", strand "<<(orfSense ? '+' : '-');
  if(headline)
    *headline=ret.str();
  ostringstream ret2;
  ret2<<'R'<<residueString<<' ';
  string shortRet;
  newCodon=origCodon;
  if(aCount) {
    newCodon[nucOffset]=orfSense ? 'A' : 'T';
    if(origCodon != newCodon) {
      ret2<<fmt2<<" A: "<<origCodon <<" -> "<<newCodon<<", "<<AminoAcidName(DNAToAminoAcid(origCodon.c_str())) <<" -> "<< AminoAcidName(DNAToAminoAcid(newCodon.c_str()))<<endl;     
    }
  }
  if(cCount) {
    newCodon[nucOffset]=orfSense ? 'C' : 'G';
    if(origCodon != newCodon) {
      ret2<<fmt2<<" C: "<<origCodon <<" -> "<<newCodon<<", "<<AminoAcidName(DNAToAminoAcid(origCodon.c_str())) <<" -> "<< AminoAcidName(DNAToAminoAcid(newCodon.c_str()))<<endl;
    }
  }
  if(gCount) {
    newCodon[nucOffset]=orfSense ? 'G' : 'C';
    if(origCodon != newCodon) {
      ret2<<fmt2<<" G: "<<origCodon <<" -> "<<newCodon<<", "<<AminoAcidName(DNAToAminoAcid(origCodon.c_str())) <<" -> "<< AminoAcidName(DNAToAminoAcid(newCodon.c_str()))<<endl;
    }
  }
  if(tCount) {
    newCodon[nucOffset]=orfSense ? 'T' : 'A';
    if(origCodon != newCodon) {
      ret2<<fmt2<<" T: "<<origCodon <<" -> "<<newCodon<<", "<<AminoAcidName(DNAToAminoAcid(origCodon.c_str())) <<" -> "<<AminoAcidName(DNAToAminoAcid(newCodon.c_str()))<<endl;
    }
  }
  if(body)
    *body=ret2.str();
  return ret2.str();
}

string makeReport(ReferenceChromosome& rg, dnapos_t pos, ReferenceChromosome::LociStats locistat, double fraction, string* summary=0)
{
  ostringstream report;
  if(summary)
    summary->clear();
  boost::format fmt1("%-10d: %3d*%c ");
  string fmt2("                  ");
  int aCount, cCount, tCount, gCount, xCount;

  vector<GeneAnnotation> gas;
  if(rg.d_gar) {
    gas= rg.d_gar->lookup("", pos);
    abort(); // needs name of chromosome here
  }

  char c=rg.snippet(pos, pos+1)[0];
  aCount = cCount = tCount = gCount = xCount = 0;
  
  acgtxDo(c, 
	 [&](){aCount += rg.d_mapping[pos].coverage;},
	 [&](){cCount += rg.d_mapping[pos].coverage;},
	 [&](){gCount += rg.d_mapping[pos].coverage;},
	  [&](){tCount += rg.d_mapping[pos].coverage;},
	  [&](){xCount += rg.d_mapping[pos].coverage;}
	 );

  char orig = rg.snippet(pos, pos+1)[0];
  report << (fmt1 % pos % rg.d_mapping[pos].coverage % orig ).str();
  sort(locistat.samples.begin(), locistat.samples.end());
  for(auto j = locistat.samples.begin(); 
      j != locistat.samples.end(); ++j) {
    c=j->nucleotide;
    acgtxDo(c, [&](){ aCount++; }, [&](){ cCount++; }, [&](){ gCount++; }, [&](){ tCount++; }, [&](){ xCount++; }  );    
    report<<c;
  }
  report<<endl<<fmt2;
  for(auto j = locistat.samples.begin(); 
      j != locistat.samples.end(); ++j) {
    report<<((char)(j->quality+33));
  }
  report << endl << fmt2;
  for(auto j = locistat.samples.begin(); 
      j != locistat.samples.end(); ++j) {
    report << (j->headOrTail ? 'H' : '.');
  }

  int tot=locistat.samples.size() + rg.d_mapping[pos].coverage;
  report<<endl;
  string aminoHeadline, aminoBody;
  if(!gas.empty())
    makeAminoReport(rg, pos, gas, locistat, &aminoHeadline, &aminoBody);
  report<<fmt2<<aminoHeadline<<endl;

  if(!gas.empty()) {
    report << fmt2 << "Annotation: ";
    for(auto& ga : gas) {
      report << ga.name<<" ["<<ga.tag<<"], ";
    }
    report << endl;
  }
  report << fmt2<< "Fraction tail: "<<fraction<<", "<< locistat.samples.size()<<endl;
  report << fmt2<< "A: " << aCount*100/tot <<"%, C: "<<cCount*100/tot<<"%, G: "<<gCount*100/tot<<"%, T: "<<tCount*100/tot<<"%"<<", X: "<<xCount*100/tot<<"%"<<endl;

  if(!aminoBody.empty())
    report << aminoBody;
  
  // cout<<rg.getMatchingFastQs(pos, fastq);
  return report.str();
}

void printQualities(FILE* jsfp, const qstats_t& qstats)
{
  int i=0;

  fprintf(jsfp, "qualities=[");
  for(const auto& q : qstats) {
    if(i)
      fputs(",", jsfp);
    if(q.valid()) {
      fprintf(jsfp, "[%d, %f]", i, -10.0*log10(mean(q)));
      ++i;
    }
  }
  fputs("];\n", jsfp);

  vector<double> qlo, qhi;
  for(const auto& q : qstats) {
    if(q.valid()) {
      qlo.push_back(-10.0*log10(mean(q)) - sqrt(-10.0*log10(variance(q))));
      qhi.push_back(-10.0*log10(mean(q)) +sqrt(-10.0*log10(variance(q))));
    }
  }

  fprintf(jsfp, "var qlo=%s;\nvar qhi=%s;\n", jsonVectorD(qlo).c_str(), jsonVectorD(qhi).c_str());

  fflush(jsfp);
}

vector<ReferenceChromosome::MatchDescriptor> getAllReadPosBoth(vector<unique_ptr<ReferenceChromosome> >& refs, const vector<unsigned int>& indexLengths, FastQRead* fqfrag) 
{
  vector<ReferenceChromosome::MatchDescriptor> ret;
  auto iter = indexLengths.begin();
  for(; iter != indexLengths.end(); ++iter)
    if(*iter == fqfrag->d_nucleotides.length())
      break;
  if(iter == indexLengths.end())
    return ret;

  for(auto& rg : refs) {
    auto inter = rg->getAllReadPosBoth(fqfrag);
    for(auto& i : inter) 
      ret.push_back(i); // XXX must be better way, back_inserter?
  }
  return ret;
}

void emitLociAndCluster(FILE* jsfp, ReferenceChromosome* rg, int numRef, 
			Clusterer<ClusterLocus>& vcl)
{
  ofstream ofs("loci."+lexical_cast<string>(numRef));
  ofs<<"locus\tnumdiff\tdepth\tA\tAq\tC\tCq\tG\tGq\tT\tTq\tdels\ttotQ\tfracHead"<<endl;
  map<dnapos_t, ReferenceChromosome::LociStats> slocimap;

  FILE* locifp=fopen(("loci."+lexical_cast<string>(numRef)+".js").c_str(), "w");
    
  for(auto& p : rg->d_locimap) {
    slocimap.insert(p);
  }


  fprintf(jsfp, "genomes[%d].loci=[", numRef);
  fprintf(locifp, "loci[\"%s\"]=[", g_name.c_str());

  bool emitted=false;
  for(auto& p : slocimap) {
    if(p.second.samples.size()==1) // no variability if only 2
      continue;
    
    int aCount=0, cCount=0, gCount=0, tCount=0, xCount=0;
    int aQual=0, cQual=0, gQual=0, tQual=0;
    int head=0;
    map<string, int> insertCounts;
    for(const auto& c: p.second.samples) {
      acgtxDo(c.nucleotide, 
	      [&](){ aCount++; aQual+=c.quality;}, 
	      [&](){ cCount++; cQual+=c.quality;}, 
	      [&](){ gCount++; gQual+=c.quality;}, 
	      [&](){ tCount++; tQual+=c.quality;},
	      [&](){ xCount++; });
      if(c.headOrTail)
	head++;
      if(!c.insert.empty())
	insertCounts[c.insert]++;
    }
    
    if(aQual < 90 && cQual < 90 && gQual < 90 && tQual < 90 && xCount < 3)
      continue;
    
    string insertReport;
    for(const auto& i : insertCounts) {
      if(!insertReport.empty())
	insertReport+=", ";
      insertReport += i.first+": "+lexical_cast<string>(i.second);
    }

    double fraction =(1.0*head/p.second.samples.size());
    if(fraction < 0.1 || fraction > 0.9)
      continue;
    
    vcl.feed(ClusterLocus{p.first, p.second});

    string summary;
    char orig = rg->snippet(p.first, p.first+1)[0];
    if(aCount && orig!='A') {
      summary.append(1, orig);
      summary.append(">A");
    }
    if(cCount && orig!='C') {
      if(!summary.empty())
	summary.append(",");
      summary.append(1, orig);
      summary.append(">C");
    }
    if(gCount && orig!='G') {
      if(!summary.empty())
	summary.append(",");
      summary.append(1, orig);
      summary.append(">G");
    }
    if(tCount && orig!='T') {
      if(!summary.empty())
	summary.append(",");
      summary.append(1, orig);
      summary.append(">T");
    }
    if(xCount) {
      if(!summary.empty())
	summary.append(",");
      summary.append("-");
      summary.append(1, orig);
    }
  

    ofs<<p.first<<"\t"<<p.second.samples.size()<<"\t"<<rg->d_mapping[p.first].coverage<<"\t";
    ofs<<aCount<<"\t"<<aQual<<"\t";
    ofs<<cCount<<"\t"<<cQual<<"\t";
    ofs<<gCount<<"\t"<<gQual<<"\t";
    ofs<<tCount<<"\t"<<tQual<<"\t";
    ofs<<xCount<<"\t";
    ofs<<aQual+cQual+gQual+tQual<<"\t";
    ofs<<fraction<<"\t"<<insertReport<<"\t#";
    string annotation;
    bool gene=false;
    string aminoReport;
    if(rg->d_gar) {
      auto gas = rg->d_gar->lookup("", p.first);
      abort(); // needs name of chromosome
      for(auto ga : gas) {
        replace_all(ga.tag, "\n", "\\n");
        replace_all(ga.tag, "'", "\\'");
	annotation+=ga.name +"\t[" + ga.tag + "]\t";
	if(ga.gene)
	  gene=1;

      }
      if(gene) {
	aminoReport = makeAminoReport(*rg, p.first, gas, p.second, 0, 0);
	replace_all(aminoReport, "\n", " ");
	trim_left(aminoReport);
      }
    }
    ofs<<annotation<<"\t"<<aminoReport<<endl;

    if(emitted) {
      fprintf(jsfp,",\n");
      fprintf(locifp, ",\n");
    }

    string graph;

    dnapos_t start = p.first-100, stop = min(p.first+100, (dnapos_t)rg->d_mapping.size());
    vector<double> aProb(stop-start), cProb(stop-start), gProb(stop-start), tProb(stop-start), xProb(stop-start);

    for(dnapos_t pos = start; pos < stop; ++pos) {
      if(!graph.empty()) 
	graph+=",";
      graph+="["+to_string(pos)+","+to_string(rg->d_mapping[pos].coverage)+"]";

      if(rg->d_locimap.count(pos)) {
	for(const auto& locus: rg->d_locimap[pos].samples) {
	  acgtxDo(locus.nucleotide, 
		  [&]() { aProb[pos-start]+= locus.quality/30.0; }, 
		  [&]() { cProb[pos-start]+= locus.quality/30.0; }, 
		  [&]() { gProb[pos-start]+= locus.quality/30.0; }, 
		  [&]() { tProb[pos-start]+= locus.quality/30.0; },
		  [&]() { xProb[pos-start]+= locus.quality/30.0; });
	}
      }
    }


    for(int n=0; n < 2; ++n) {
      fprintf(n ? jsfp : locifp, " { locus: %d, numDiff: %d, originalBase: '%c', depth: %d, "
	      "aCount: %d, aQual: %d, "
	      "cCount: %d, cQual: %d, "
	      "gCount: %d, gQual: %d, "
	      "tCount: %d, tQual: %d, "
	      "totQual: %d, "
	      "xCount: %d, "
	      "fraction: %f, gene: %d, annotation: '%s', aminoReport: '%s', insertReport: '%s', summary: '%s', graph: [%s], aProb: %s, cProb: %s, gProb: %s, tProb: %s, xProb: %s}", 
	      p.first, (int)p.second.samples.size(), '?', rg->d_mapping[p.first].coverage, 
	      aCount, aQual, cCount, cQual, gCount, gQual, tCount, tQual, 
	      aQual+cQual+gQual+tQual,
	      xCount, fraction, gene, annotation.c_str(), aminoReport.c_str(), insertReport.c_str(), summary.c_str(), graph.c_str(),
	      jsonVectorX(aProb, [start](int i){return i+start;}).c_str(), 
	      jsonVectorX(cProb, [start](int i){return i+start;}).c_str(), 
	      jsonVectorX(gProb, [start](int i){return i+start;}).c_str(), 
	      jsonVectorX(tProb, [start](int i){return i+start;}).c_str(),
	      jsonVectorX(xProb, [start](int i){return i+start;}).c_str());
    }
    emitted=true;
  }
  fprintf(jsfp, "];\n");
  fprintf(locifp, "];\n");
  fclose(locifp);
  ofs.flush();
}

template<typename T>
void safeIncVec(vector<T>& vec, unsigned int offset) 
{
  if(offset >= vec.size())
    vec.resize(offset+1);
  vec[offset]++;
}

bool g_pleaseQuit;
void pleaseQuitHandler(int)
{
  g_pleaseQuit=true;
}


void doInitialReadStatistics(FILE* jsfp, const string& fname, StereoFASTQReader& fastq, vector<unsigned int>* recommendIndex, unsigned int* maxreadlen, unsigned int *recommendBeginSnip=0, unsigned int* recommendEndSnip=0)
{
  FastQRead fqfrag1, fqfrag2;
  vector<uint32_t> lengths;

  vector<vector<uint32_t>> kmerMappings(1024);
  vector<uint32_t> gcMappings(1024), taMappings(1024);
  
  for(auto& kmers : kmerMappings) 
    kmers.resize(256); // 4^4, corresponds to the '4' below
  
  (*g_log)<<"Scanning FASTQ input to determine trim optima and indexation parameters"<<endl;
  uint64_t totalReads=0;

  boost::progress_display show_progress(filesize(fname.c_str()), cerr);
  unsigned int bytes;
  unsigned int counter=0;
  while((bytes=fastq.getReadPair(&fqfrag1, &fqfrag2))) {
    show_progress+=bytes;
    if((++counter%11)) 
      continue;
    safeIncVec(lengths, fqfrag1.d_nucleotides.length());
    safeIncVec(lengths, fqfrag2.d_nucleotides.length());
    totalReads+=2;

    for(int n=0; n < 2; ++n) {
      FastQRead* fqfrag = n ? &fqfrag1 : &fqfrag2;
      for(string::size_type i = 0 ; i < fqfrag->d_nucleotides.size(); ++i) {
	if(fqfrag->d_nucleotides.size() - i > 4)
	  kmerMappings[i][kmerMapper(fqfrag->d_nucleotides, i, 4)]++;

	char c = fqfrag->d_nucleotides[i];
	if(c=='G' || c=='C')
	  gcMappings[i]++;  
	else
	  taMappings[i]++;
      }
    }
  }

  multimap<uint32_t, uint32_t> lenstats;
  *maxreadlen=0;
  for(auto iter = lengths.cbegin(); iter != lengths.cend(); ++iter) {
    lenstats.insert(make_pair(*iter, iter - lengths.cbegin())); // amount, length
    *maxreadlen = max(*maxreadlen, (unsigned int)(iter - lengths.cbegin()));
  }

  uint64_t cumul=0;
  int num=0;
  for(auto iter = lenstats.crbegin(); iter != lenstats.crend() && num < 3; ++iter, ++num) {
    cumul+=iter->first;
    (*g_log) << cumul*100.0/totalReads<<"% covered after indexing "<<iter->second<<endl;
    recommendIndex->push_back(iter->second);
    if(cumul > 0.99*totalReads)
      break;
  }
  fastq.seek(0);

  fprintf(jsfp, "var kmerstats=[");
  unsigned int readOffset=0;
  for(const auto& kmer : kmerMappings) {
    if(readOffset >= kmerMappings.size() - 4)
      break;
    
    VarMeanEstimator acc;
    for(auto& count : kmer) {
      acc(count);
    }
    if(mean(acc)!=0)
      fprintf(jsfp, "%s[%d, %f]", readOffset ? "," : "", readOffset, sqrt(variance(acc)) / mean(acc) );
    readOffset++;
  }
  fprintf(jsfp, "];\n");
  
  vector<double> ratios;
  VarMeanEstimator ratest;
  fprintf(jsfp, "var gcRatios=[");
  for(unsigned int i=0; i < *maxreadlen ;++i) {
    double total=0.001+gcMappings[i] + taMappings[i];
    double ratio= gcMappings[i]/total;

    fprintf(jsfp, "%s[%d,%.2f]", i ? "," : "", i, ratio);
    
    ratios.push_back(ratio);
    if(i > 0.1 * *maxreadlen && i < 0.9 * *maxreadlen)
      ratest(ratio);
  }
  fprintf(jsfp,"];\n");

  fflush(jsfp);

  //  cout<<"Variance: "<<sqrt(variance(ratest))<<endl;
  // so where do we put the cut..
  for(unsigned int n=0; n < ratios.size(); ++n) {
    //    cout<<n<<" "<< (ratios[n]-mean(ratest))/sqrt(variance(ratest))<<endl;
    int j=10;
    for(; j  && n+j < ratios.size(); --j) {
      double sigma=fabs((ratios[n+j]-mean(ratest))/sqrt(variance(ratest)));
      if(sigma > 5.0)
	break;
      //      cout<<" "<<n+j<<" "<< (ratios[n+j]-mean(ratest))/sqrt(variance(ratest))<<endl;
    }
    if(!j) {
      (*g_log)<<"Put the begin trim at: "<<n+1<<endl;
      if(recommendBeginSnip) {
	*recommendBeginSnip=n+1;
	for(auto& i: *recommendIndex)
	  i-=*recommendBeginSnip;
	*maxreadlen-=*recommendBeginSnip;
      }
      break;
    }
  }
}


int main(int argc, char** argv)
try
{
#ifdef __linux__
  feenableexcept(FE_DIVBYZERO | FE_INVALID); 
#endif 
  srand(time(0));  
  TCLAP::CmdLine cmd("Command description message", ' ', "g" + string(g_gitHash));

  TCLAP::MultiArg<std::string> annotationsArg("a","annotations","read annotations for reference genome from this file",false, "filename", cmd);
  TCLAP::MultiArg<std::string> referenceArg("r","reference","read annotations for reference genome from this file",true,"string", cmd);
  TCLAP::ValueArg<std::string> fastq1Arg("1","fastq1","read annotations for reference genome from this file",true,"","string", cmd);
  TCLAP::ValueArg<std::string> fastq2Arg("2","fastq2","read annotations for reference genome from this file",true,"","string", cmd);

  TCLAP::ValueArg<std::string> nameArg("n","name","Name for this analysis",false,"","string", cmd);

  TCLAP::ValueArg<std::string> bamFileArg("w","bam-file","Write the assembly to the named BAM file",false,"","filename", cmd);
  TCLAP::ValueArg<int> qualityOffsetArg("q","quality-offset","Quality offset in fastq. 33 for Sanger.",false, 33,"offset", cmd);
  TCLAP::ValueArg<int> beginSnipArg("b","begin-snip","Number of nucleotides to snip from begin of reads",false, 0,"nucleotides", cmd);
  TCLAP::ValueArg<int> endSnipArg("e","end-snip","Number of nucleotides to snip from end of reads",false, 0,"nucleotides", cmd);
  TCLAP::ValueArg<int> qlimitArg("l","qlimit","Disregard nucleotide reads with less quality than this in calls",false, 30,"q", cmd); 
  TCLAP::ValueArg<int> duplimitArg("d","duplimit","Ignore reads that occur more than d times. 0 for no filter.",false, -1,"times", cmd);
  TCLAP::SwitchArg unmatchedDumpSwitch("u","unmatched-dump","Create a dump of unmatched reads (unfound.fastq)", cmd, false);
  TCLAP::SwitchArg skipUndermatchedSwitch("","skip-undermatched","Do not emit undermatched regions", cmd, true);
  TCLAP::SwitchArg skipVariableSwitch("","skip-variable","Do not emit variable regions", cmd, false);
  TCLAP::SwitchArg skipInsertsSwitch("","skip-inserts","Do not emit inserts", cmd, false);
  TCLAP::SwitchArg excludePhiXSwitch("p","exclude-phix","Exclude PhiX automatically",cmd, false);

  cmd.parse( argc, argv );

  unsigned int qlimit = qlimitArg.getValue();

  ostringstream jsonlog;  
  TeeDevice td(cerr, jsonlog);
  g_log = new TeeStream(td);

  (*g_log)<<"Antonie was compiled from git hash g" << g_gitHash <<", using "<<compilerVersion() << " on " << __TIME__ << " " <<__DATE__ <<endl;
  g_name=nameArg.getValue();
  if(!g_name.empty()) {
    (*g_log)<<"Name of this run: "<<g_name<<endl;  
  }

  {
    char buffer[1024];
    if(getcwd(buffer,sizeof(buffer)))
      (*g_log)<<"Current working directory is '"<<buffer<<"' on host ";
    *buffer=0;
#ifdef _WIN32
    DWORD len=sizeof(buffer);
    if(GetComputerName(buffer, &len)) {
      buffer[len]=0;
      (*g_log)<<buffer << endl;
    }
    else (*g_log)<<"UNKNOWN"<<endl;
#else 
    (*g_log)<<(gethostname(buffer, sizeof(buffer)) < 0 ? "UNKNOWN" : buffer) <<endl;
#endif
  }
  //  (*g_log)<<"Current time: "<< std::put_time(std::localtime(&system_clock::now()), "%F %T")<<endl;
  
  StereoFASTQReader fastq(fastq1Arg.getValue(), fastq2Arg.getValue(), qualityOffsetArg.getValue());

  (*g_log)<<"FASTQ Input from '"<<fastq1Arg.getValue()<<"' and '"<<fastq2Arg.getValue()<<"'"<<endl;
  unique_ptr<FILE, int(*)(FILE*)> jsfp(fopen("data.js","w"), fclose);
  vector<unsigned int> indexLengths;
  unsigned int maxreadsize=0;
  unsigned int beginTrim=beginSnipArg.getValue(), endTrim= endSnipArg.getValue();
  doInitialReadStatistics(jsfp.get(), fastq1Arg.getValue(), fastq, &indexLengths, &maxreadsize, beginTrim ? 0 :&beginTrim, endTrim ? 0 : &endTrim);
  fastq.setTrim(beginTrim, endTrim);
  (*g_log)<<"Trimming "<<beginTrim<<" from beginning of reads, "<<endTrim<<" from end of reads"<<endl;

  unsigned int bytes=0;
  FastQRead fqfrag1, fqfrag2;


  bytes=fastq.getReadPair(&fqfrag1, &fqfrag2); // get a read to index based on its size
  int keylen=11;

  fputs("var genomes=[];\n", jsfp.get());

  vector<unique_ptr<ReferenceChromosome> > refgens;
  auto annotations = annotationsArg.getValue().begin();
  for(auto& fname: referenceArg.getValue()) {
    unique_ptr<ReferenceChromosome> rg(new ReferenceChromosome(fname));
    double genomeGCRatio = 1.0*(rg->d_cCount + rg->d_gCount)/(rg->d_cCount + rg->d_gCount + rg->d_aCount + rg->d_tCount);

    (*g_log)<<"Read FASTA reference genome of '"<<rg->d_fullname<<"', "<<rg->size()<<" nucleotides from '"<<fname<<"' (GC = "<<genomeGCRatio<<")"<<endl;
    for(auto i : indexLengths)
      rg->index(i);
    rg->index(keylen);
    fprintf(jsfp.get(), "var genomeGCRatio=%f;\n", genomeGCRatio); // XXXmulti

    if(annotations != annotationsArg.getValue().end()) {
      auto gar = new GeneAnnotationReader(*annotations);
      (*g_log)<<"Done reading "<<gar->size()<<" annotations from '"<<*annotations<<"'"<<endl;
      rg->addAnnotations(gar);  
      annotations++;
    }
    else
      (*g_log)<<"No annotations for '"<<rg->d_fullname<<"': "<<((bool)rg->d_gar)<<endl;
    refgens.emplace_back(move(rg));
  }

  if(excludePhiXSwitch.getValue()) {
    auto rg = ReferenceChromosome::makeFromString(phiXFastA);
    double genomeGCRatio = 1.0*(rg->d_cCount + rg->d_gCount)/(rg->d_cCount + rg->d_gCount + rg->d_aCount + rg->d_tCount);
    (*g_log)<<"Read FASTA reference genome of '"<<rg->d_fullname<<"', "<<rg->size()<<" nucleotides from builtin (GC = "<<genomeGCRatio<<")"<<endl;
    rg->index(maxreadsize);
    rg->index(keylen);

    auto gar = new GeneAnnotationReader("./phix.gff");
    (*g_log)<<"Done reading "<<gar->size()<<" annotations from builtin"<<endl;
    rg->addAnnotations(gar);  

    refgens.emplace_back(move(rg));
  }

  int duplimit = duplimitArg.getValue();
  if(duplimit < 0) {
    dnapos_t totsize=0;
    for(auto& rg : refgens) 
      totsize += rg->size();
    auto lambda = 1.0*fastq.estimateReads()/totsize;
    boost::math::poisson_distribution<double> pd(lambda);
    duplimit=boost::math::quantile(pd, 0.999);
    (*g_log)<<"Auto-set duplicate filter to "<<duplimit<<" based on 0.999 cumulative Poisson. Expect "<<lambda<<" dups on average over estimated "<< fastq.estimateReads()<<" reads"<<endl;
  }
  else if(duplimit > 0)
    (*g_log)<<"Duplicate reads filtered beyond "<<duplimit<<" copies"<<endl;

  g_log->flush();
  dnapos_t pos;

  uint64_t withAny=0, found=0, total=0, tooFrequent=0, goodPairMatches=0, badPairMatches=0,
    qualityExcluded=0;

  BAMWriter sbw(bamFileArg.getValue(), (*refgens.begin())->d_name, (*refgens.begin())->size()); // XXXmulti

  (*g_log)<<"Performing matches of reads to reference genome"<<endl;
  boost::progress_display show_progress(filesize(fastq1Arg.getValue().c_str()), cerr);
 

  qstats_t qstats;
  qstats.resize(maxreadsize);
  VarMeanEstimator qstat;
  vector<unsigned int> qcounts(256);
  vector<uint64_t> unfoundReads;
  vector<qtally> qqcounts(256);
  vector<dnapos_t> gchisto(maxreadsize+1);

  DuplicateCounter dc;
  uint32_t theHash;
  map<uint32_t, uint32_t> seenAlready;
  vector<uint32_t> pairdisthisto;
  vector<uint32_t> readlengths;
  signal(SIGINT, pleaseQuitHandler);

  do { 
    if(g_pleaseQuit) 
      break;
    show_progress += bytes;
    vector<ReferenceChromosome::MatchDescriptor > pairpositions[2];
    bool dup1(false), dup2(false);
    safeIncVec(readlengths, fqfrag1.d_nucleotides.length());
    safeIncVec(readlengths, fqfrag2.d_nucleotides.length());
    for(unsigned int paircount=0; paircount < 2; ++paircount) {
      FastQRead& fqfrag(paircount ? fqfrag2 : fqfrag1);
      total++;
      for(string::size_type pos = 0 ; pos < fqfrag.d_quality.size(); ++pos) {
	int i = fqfrag.d_quality[pos];
	double err = qToErr(i);
	qstat(err);
	qstats[pos](err);
	qcounts[i]++;
      }
      dc.feedString(fqfrag.d_nucleotides);
      if(duplimit) {
	theHash=qhash(fqfrag.d_nucleotides.c_str(), fqfrag.d_nucleotides.size(), 0);
	if(++seenAlready[theHash] > (unsigned int)duplimit) {
	  if(paircount)
	    dup2=true;
	  else
	    dup1=true;
	  tooFrequent++;
	  continue;
	}
      }
      
      gchisto[round(fqfrag.d_nucleotides.size()*getGCContent(fqfrag.d_nucleotides))]++;
      
      if(fqfrag.d_nucleotides.find('N') != string::npos) {
	// unfoundReads.push_back(fqfrag.position); // will fail elsewhere and get filed there
	withAny++;
	continue;
      }
      if((pairpositions[paircount]=getAllReadPosBoth(refgens, indexLengths, &fqfrag)).empty()) {
	pairpositions[paircount]=fuzzyFind(&fqfrag, refgens, keylen, qlimit);
      } 
    }
    
    map<int, vector<pair<ReferenceChromosome::MatchDescriptor,ReferenceChromosome::MatchDescriptor> > > potMatch;
    unsigned int matchCount=0;
    for(auto& match1 : pairpositions[0]) {
      for(auto& match2 : pairpositions[1]) {
	if(match1.reverse != match2.reverse &&  abs((int64_t) match1.pos - (int64_t)match2.pos) < 1400) {
	  potMatch[match1.score + match2.score].push_back({match1, match2});
	  matchCount++;
	}
      }
    }
   
    if(!potMatch.empty()) {
      const auto& chosen = pickRandom(potMatch.begin()->second);
      goodPairMatches++;
      int distance = chosen.second.reverse ? 
	(fqfrag1.d_nucleotides.length() + (int64_t) chosen.second.pos - (int64_t) chosen.first.pos) :
	(fqfrag1.d_nucleotides.length() + (int64_t) chosen.first.pos - (int64_t) chosen.second.pos);

      if(distance >= 0)
	safeIncVec(pairdisthisto, distance);
      for(int paircount = 0 ; paircount < 2; ++paircount) {
	auto fqfrag = paircount ? &fqfrag2 : &fqfrag1;
	auto dup = paircount ? dup2 : dup1,
	  otherDup = paircount? dup1 : dup2;
	pos = paircount ? chosen.second.pos : chosen.first.pos;


	if((paircount ? chosen.second.reverse : chosen.first.reverse) != fqfrag->reversed)
	  fqfrag->reverse();

	if(otherDup && !dup) {
	  MapToReference(*chosen.second.rg, pos, *fqfrag, qlimit, &sbw, &qqcounts);
	}
	else if(!otherDup && !dup) {
	  int indel; 
	  // XXX add amount here
	  if(MapToReference(*chosen.second.rg, pos, *fqfrag, qlimit, 0, &qqcounts, &indel)) {
	    sbw.qwrite(pos, *fqfrag, indel, 3 + (paircount ? 0x80 : 0x40),
		     "=", 
		     paircount ? chosen.first.pos : chosen.second.pos, 
		       (chosen.first.reverse ^ (bool)paircount) ? -distance : distance);
	  }
	}
	found++;
      }
    }
    else {
      //      cout<<"No pair matches, need to map individually: "<<endl;
      badPairMatches++;
      for(unsigned int paircount = 0; paircount < 2; ++paircount) {
	if(paircount ? dup2 : dup1)
	  continue;
	
	map<int, vector<ReferenceChromosome::MatchDescriptor>> scores;
	for(auto match: pairpositions[paircount]) {
	  scores[match.score].push_back(match);
	  //	  cout<<"\t"<<paircount<<"\t"<<match.pos<<" "<<match.reverse<<", score: "<<match.score<<endl;
	}
	FastQRead* fqfrag = paircount ? &fqfrag2 : &fqfrag1;
	if(scores.empty()) {
	  unfoundReads.push_back(fqfrag->position);
	  continue;
	}
	auto pick = pickRandom(scores.begin()->second);

	if(fqfrag->reversed != pick.reverse)
	  fqfrag->reverse();

	MapToReference(*pick.rg, pick.pos, *fqfrag, qlimit, &sbw, &qqcounts);
	found++;
      }
    } 
  } while((bytes=fastq.getReadPair(&fqfrag1, &fqfrag2)));
  signal(SIGINT, SIG_DFL);
  
  pairdisthisto.resize(1500); // outliers mess us up otherwise
  fputs(jsonVector(pairdisthisto, "var pairdisthisto").c_str(), jsfp.get());
  fputs(jsonVector(readlengths, "var readlengths").c_str(), jsfp.get());

  uint64_t totNucleotides=total*maxreadsize; // XXX very wrong
  fprintf(jsfp.get(), "qhisto=[");
  for(int c=0; c < 50; ++c) {
    fprintf(jsfp.get(), "%s[%d,%f]", c ? "," : "", (int)c, 1.0*qcounts[c]/totNucleotides);
  }
  fprintf(jsfp.get(),"];\n");

  fprintf(jsfp.get(), "var dupcounts=[");
  auto duplicates = dc.getCounts();
  for(auto iter = duplicates.begin(); iter != duplicates.end(); ++iter) {
    fprintf(jsfp.get(), "%s[%" PRIu64 ",%f]", (iter!=duplicates.begin()) ? "," : "", iter->first, 1.0*iter->second/total);
  }
  fprintf(jsfp.get(),"];\n");
  dc.clear(); // might save some memory..

  dnapos_t totalhisto= accumulate(gchisto.begin(), gchisto.end(), 0);
  fputs(jsonVector(gchisto, "var gcreadhisto",  
		   [totalhisto](dnapos_t c){return 1.0*c/totalhisto;},
		   [&maxreadsize](int i) { return 100.0*i/maxreadsize;}   ).c_str(),  // XXX wrong scaling
	jsfp.get());

  unsigned int numRef=0;
  for(auto& rg : refgens) {
    fprintf(jsfp.get(), "genomes[%d]={};\n", numRef);
    fputs(jsonVector(rg->getGCHisto(), ("genomes["+lexical_cast<string>(numRef)+"].gcrefhisto").c_str(),  
		     [&maxreadsize,&rg](dnapos_t c){return 1.0*c/(rg->size()/maxreadsize);},
		     [&maxreadsize](int i) { return 100.0*i/maxreadsize;}   ).c_str(),  // XXX wrong scaling
	  jsfp.get());


    numRef++;
  }

  (*g_log) << (boost::format("Total reads: %|40t| %10d (%.2f gigabps)") % total % (totNucleotides/1000000000.0)).str() <<endl;
  (*g_log) << (boost::format("Quality excluded: %|40t|-%10d") % qualityExcluded).str() <<endl;
  (*g_log) << (boost::format("Ignored reads with N: %|40t|-%10d") % withAny).str()<<endl;
  if(duplimit)
    (*g_log) << (boost::format("Too frequent reads: %|40t| %10d (%.02f%%)") % tooFrequent % (100.0*tooFrequent/total)).str() <<endl;
  (*g_log) << (boost::format("Full matches: %|40t|-%10d (%.02f%%)\n") % found % (100.0*found/total)).str();
  (*g_log) << (boost::format(" Reads matched in a good pair: %|40t| %10d\n") % (goodPairMatches*2)).str();
  (*g_log) << (boost::format(" Reads not matched, bad pair: %|40t| %10d\n") % (badPairMatches*2)).str();

  (*g_log) << (boost::format("Not fully matched: %|40t|=%10d (%.02f%%)\n") % unfoundReads.size() % (unfoundReads.size()*100.0/total)).str();
  (*g_log) << (boost::format("Mean Q: %|40t|    %10.2f +- %.2f\n") % (-10.0*log10(mean(qstat))) 
	       % sqrt(-10.0*log10(variance(qstat)) )).str();

  seenAlready.clear();

  for(auto& rg : refgens) {  // XXXmulti - the 'found' should be per GC, not global!
    for(auto& i : rg->d_correctMappings) {
      i=found;
    }
  }
  printQualities(jsfp.get(), qstats);

  if(!bamFileArg.getValue().empty()) {
    (*g_log) << "Writing sorted & indexed BAM file to '"<< bamFileArg.getValue()<<"'"<<endl;
    sbw.runQueue(fastq);
  }
  if(unmatchedDumpSwitch.getValue())
    writeUnmatchedReads(unfoundReads, fastq);
  int index=0;
  numRef=0;

  for(auto& rg : refgens) {
    (*g_log)<<"Output for "<<rg->d_fullname<<endl;
    rg->printCoverage(jsfp.get(), "genomes["+lexical_cast<string>(numRef)+"].fullHisto");
    Clusterer<Unmatched> cl(100);
    for(auto unm : rg->d_unmRegions) {
      cl.feed(unm);
    }
    if(skipUndermatchedSwitch.getValue()) {
      (*g_log)<<"Skipping output of undermatched regions"<<endl;
    }
    else {
      for(auto unmCl : cl.d_clusters) {
	string report=makeReport(*rg, pos, rg->d_locimap[pos], -1);
	emitRegion(jsfp.get(), *rg, fastq, "Undermatched", index++, unmCl.getBegin()-100, unmCl.getEnd()+100, report);
      }
    }
    printCorrectMappings(jsfp.get(), *rg, "genomes["+lexical_cast<string>(numRef)+"].referenceQ");
    fprintf(jsfp.get(), "genomes[%d].qqdata=[", numRef);
    bool printedYet=false;
    for(auto coinco = qqcounts.begin() ; coinco != qqcounts.end(); ++coinco) {
      if(coinco->incorrect || coinco->correct) {
	double qscore;
	if(coinco->incorrect && coinco->correct)
	  qscore = -10.0*log10(1.0*coinco->incorrect / (coinco->correct + coinco->incorrect));
	else if(coinco->correct == 0)
	  qscore=0;
	else
	  qscore=41; // "highest score possible"
	
	fprintf(jsfp.get(), "%s[%u, %f, %" PRIu64 "]", 
		printedYet ?  "," : "", 
		(unsigned int)(coinco - qqcounts.begin()), qscore,
		coinco->incorrect + coinco->correct);
	printedYet=true;
      }
    }
    fprintf(jsfp.get(), "];\n");

    struct revsort
    {
      bool operator()(const unsigned int&a, const unsigned int&b) const
      { return a > b;} 
    };
    (*g_log)<<"Found "<<rg->d_insertCounts.size()<<" loci with at least one insert in a read"<<endl;
    map<unsigned int, vector<dnapos_t>, revsort> topInserts;
    unsigned int significantInserts=0;
    for(const auto& insloc : rg->d_insertCounts) {
      topInserts[insloc.second].push_back(insloc.first);
      if(insloc.second > 4)
	significantInserts++;
    }
    (*g_log)<<"Found "<<significantInserts<<" significant inserts"<<endl;


    uint64_t significantlyVariable=0;
    Clusterer<ClusterLocus> vcl(100);
    emitLociAndCluster(jsfp.get(), rg.get(), numRef, vcl);
    (*g_log)<<vcl.numClusters()<<" clusters of real variability, " << vcl.numEntries()<<" variable loci"<<endl;
    
    if(skipVariableSwitch.getValue()) {
      (*g_log)<<"Not emitting variable regions"<<endl;
    }
    else {
      for(auto& cluster : vcl.d_clusters) {
	vector<pair<string, dnapos_t>> reports;
	int maxVarcount=0;
	for(auto& locus : cluster.d_members) {
	  double fraction=0;
	  int varcount=variabilityCount(*rg, locus.pos, locus.locistat, &fraction);
	  maxVarcount = max(varcount, maxVarcount);
	  //if(varcount < 3) 
	  //  continue;
	  significantlyVariable++;
	  string report = makeReport(*rg, locus.pos, locus.locistat, fraction);
	  reports.push_back({report, locus.pos});  
	}
	if(!reports.empty()) {
	  string theReport;
	  dnapos_t minPos=reports.begin()->second;
	  dnapos_t maxPos=minPos;
	  for(auto r : reports) {
	    theReport += r.first;
	    minPos=min(r.second, minPos);
	    maxPos=max(r.second, maxPos);
	  }	  

	  emitRegion(jsfp.get(), *rg, fastq, "Variable", index++, minPos > 100 ? minPos-100 : 1, maxPos+100 > rg->size() ? rg->size() : maxPos+100, theReport, maxVarcount);
	}
      }
      (*g_log)<<"Found "<<significantlyVariable<<" significantly variable loci"<<endl;
    }
    if(skipInsertsSwitch.getValue()) {
      (*g_log)<<"Not emitting inserts"<<endl;
    }
    else {
      for(const auto& insert : topInserts) {
	if(insert.first < 3)
	  break;
	for(const auto& position : insert.second) {
	  auto theReport = makeReport(*rg, position, rg->d_locimap[position], 0);
	  emitRegion(jsfp.get(), *rg, fastq, "Insert", index++, position, theReport);
	}
      }
    }
    numRef++;
  }
  g_log->flush();
  string log = jsonlog.str();
  replace_all(log, "\n", "\\n");
  fprintf(jsfp.get(), "var antonieLog=\"%s\";\n", log.c_str());

  return EXIT_SUCCESS;
}
catch(exception& e)
{
  cerr<<"Had error: "<<e.what()<<endl;
}
