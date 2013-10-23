#define __STDC_FORMAT_MACROS
#include <stdio.h>
#include <iostream>
#include <stdint.h>
#include <map>
#include <vector>
#include <string>
#include <string.h>
#include <stdexcept>
#include <inttypes.h>
#include <boost/progress.hpp>
#include <algorithm>
#include <numeric>
#include <unistd.h>
#include <errno.h>
#include <math.h>
#include <boost/format.hpp>
#include <array>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/tuple/tuple.hpp>
#include <memory>
#include "geneannotated.hh"
#include "misc.hh"
#include "fastq.hh"
// #include <thread>
#include <mba/diff.h>
#include <mba/msgno.h>

using namespace std;
using namespace boost::accumulators;

extern "C" {
#include "hash.h"
}

typedef uint32_t dnapos_t;
unsigned int dnanpos = (uint32_t) -1;

struct FASTQMapping
{
  uint64_t pos;
  bool reverse;
  int indel; // 0 = nothing, >0 means WE have an insert versus reference at pos
             // <0 means WE have a delete versus reference at pos
};

struct GenomeLocusMapping
{
  GenomeLocusMapping() : coverage(0) {}
  vector<FASTQMapping> d_fastqs;
  unsigned int coverage;
};
 
class ReferenceGenome
{
public:
  ReferenceGenome(const string& fname);
  dnapos_t size() const {
    return d_genome.size();
  }
  vector<uint32_t> getMatchingHashes(const vector<uint32_t>& hashes);
  vector<dnapos_t> getReadPositions(const std::string& nucleotides)
  {
    vector<dnapos_t> ret;
    if(nucleotides.length() != d_indexlength)
      throw runtime_error("Attempting to find a read of length we've not indexed for");
      
    uint32_t hashval = hash(nucleotides.c_str(), nucleotides.length(), 0);
    HashPos hp(hashval, 0);
    pair<index_t::const_iterator, index_t::const_iterator> range = equal_range(d_index.begin(), d_index.end(), hp);
    if(range.first == range.second)
      return ret;

    for(;range.first != range.second; range.first++) {
      if(!memcmp(d_genome.c_str() + range.first->d_pos, nucleotides.c_str(), nucleotides.length())) {
        ret.push_back(range.first->d_pos);
      }
    }
    return ret;
  }
  
  dnapos_t getReadPosBoth(FastQRead* fq) // tries both
  {
    vector<dnapos_t> positions;
    string nucleotides;
    for(int tries = 0; tries < 2; ++tries) {
      positions = getReadPositions(fq->d_nucleotides);

      if(!positions.empty()) {
        int rpos=random() % positions.size();
        cover(positions[rpos], fq->d_nucleotides.size(), fq->d_quality);
        return positions[rpos];
      }
      fq->reverse();
    }
    return dnanpos;
  }

  void cover(dnapos_t pos, unsigned int length, const std::string& quality) 
  {
    const char* p = quality.c_str();
    for(unsigned int i = 0; i < length; ++i) {
      if(p[i]-33 > 30)
        d_mapping[pos+i].coverage++;
    }
  }

  void cover(dnapos_t pos, char quality) 
  {
    if(quality-33 > 30)
      d_mapping[pos].coverage++;
  }

  void mapFastQ(dnapos_t pos, const FastQRead& fqfrag, int indel=0)
  {
    FASTQMapping fqm;
    fqm.pos=fqfrag.position;
    fqm.reverse = fqfrag.reversed;
    fqm.indel = indel;
    d_mapping[pos].d_fastqs.push_back(fqm);
    //    cout<<"Adding mapping at pos "<<pos<<", indel = "<<indel<<", reverse= "<<fqm.reverse<<endl;
  }

  string snippet(dnapos_t start, dnapos_t stop) const { 
    if(stop > d_genome.size()) {
      return d_genome.substr(start);
    }
    return d_genome.substr(start, stop-start);
  }
  void printCoverage(FILE* jsfp);
  void index(int length);
  void printFastQs(dnapos_t pos, FASTQReader& fastq);
  vector<GenomeLocusMapping> d_mapping;
  vector<unsigned int> d_correctMappings, d_wrongMappings;
private:
  unsigned int d_indexlength;
  string d_genome;


  struct HashPos {
    HashPos(uint32_t hash_, dnapos_t pos) : d_hash(hash_), d_pos(pos)
    {}
    HashPos(){}
    uint32_t d_hash;
    dnapos_t d_pos;
    
    bool operator<(const HashPos& rhs) const 
    {
      return d_hash < rhs.d_hash;
    }
  };

  typedef vector<HashPos> index_t;
  index_t d_index;
};

ReferenceGenome::ReferenceGenome(const string& fname)
{
  FILE* fp = fopen(fname.c_str(), "r");
  if(!fp)
    throw runtime_error("Unable to open reference genome file '"+fname+"'");
  d_genome.reserve(filesize(fname.c_str()));  // slight overestimate which is great
  d_correctMappings.resize(150);
  d_wrongMappings.resize(150);

  char line[256]="";

  sfgets(line, sizeof(line), fp);
  chomp(line);

  if(line[0] != '>') 
    throw runtime_error("Input not FASTA");
  cerr<<"Reading FASTA reference genome of '"<<line+1<<"'\n";

  d_genome="*"; // this gets all our offsets ""right""
  while(fgets(line, sizeof(line), fp)) {
    chomp(line);
    d_genome.append(line);
  }
  
  d_mapping.resize(d_genome.size());
  d_indexlength=0;
}

void ReferenceGenome::index(int length)
{
  d_index.clear();
  d_index.reserve(d_genome.length());
  d_indexlength=length;
  cerr<<"Indexing "<<d_genome.length()<<" nucleotides for length "<<length<<"..";

  //  boost::progress_display show_progress( d_genome.length() - d_indexlength, cerr);

  for(string::size_type pos = 0 ; pos < d_genome.length() - length; ++pos) {
    uint32_t hashval = hash(d_genome.c_str() + pos, d_indexlength, 0);
    d_index.push_back(HashPos(hashval, pos));
    // ++show_progress;
  }

  cerr<<" done\nSorting hashes..";
  sort(d_index.begin(), d_index.end());
  cerr<<" done"<<endl;
    
  //  cerr<<"Average hash fill: "<<1.0*d_genome.length()/d_index.size()<<endl;
}

void ReferenceGenome::printFastQs(dnapos_t pos, FASTQReader& fastq) 
{
  if(pos < 150)
    pos = 0;
  else
    pos -= 150;
  
  string reference=snippet(pos, pos+300);
  unsigned int insertPos=0;
  for(unsigned int i = 0 ; i < 251; ++i) {
    if(i==75)
      cout << reference << endl;
    string spacer(i, ' ');
    for(auto& fqm : d_mapping[pos+i].d_fastqs) {
      fastq.seek(fqm.pos);
      FastQRead fqr;
      fastq.getRead(&fqr);
      if(fqm.reverse)
        fqr.reverse();

      if(fqm.indel > 0 && !insertPos) { // our read has an insert at this position, stretch reference
        if(i+fqm.indel < reference.size())
          reference.insert(i+fqm.indel, 1, '_');
        insertPos=i+fqm.indel;
      } else if(fqm.indel < 0) {      // our read has an erase at this position
        fqr.d_nucleotides.insert(-fqm.indel, 1, 'X');
        fqr.d_quality.insert(-fqm.indel, 1, 'X');
      }
      
      if(fqm.indel <= 0 && insertPos && i > insertPos) {
        fqr.d_nucleotides.insert(0, 1, '<');
        fqr.d_quality.insert(0, 1, 'X');
      }
      cout << spacer;
      int offset=0;
      for(unsigned int j = 0 ; j < fqr.d_nucleotides.size() && i + j + offset < reference.size(); ++j) {
        if(reference[i+j]=='_' && !fqm.indel) {
          cout<<'_';
          offset=1;
        }
        if(reference[i+j+offset]==fqr.d_nucleotides[j])
          cout<<'.';
        else if(fqr.d_quality[j] > '@') 
          cout << fqr.d_nucleotides[j];
        else if(fqr.d_quality[j] < '7') 
          cout << ' ';
        else
          cout<< (char)tolower(fqr.d_nucleotides[j]);
      }
      cout<<"                 "<<(fqm.reverse ? 'R' : ' ');
      cout<<endl;
    }
  }
}

struct Unmatched
{
  string left, unmatched, right;
  dnapos_t pos;
};


vector<uint32_t> ReferenceGenome::getMatchingHashes(const vector<uint32_t>& hashes)
{
  struct cmp
  {
    bool operator()(const uint32_t& a, const HashPos& b) const
    {
      return a < b.d_hash;
    }
    bool operator()(const HashPos& a, const uint32_t& b) const
    {
      return a.d_hash < b;
    }

  };
  vector<uint32_t> intersec;
  set_intersection(hashes.begin(), hashes.end(), d_index.begin(), d_index.end(), back_inserter(intersec), cmp());

  return intersec;
}
vector<Unmatched> g_unm;

void ReferenceGenome::printCoverage(FILE* jsfp)
{
  uint64_t totCoverage=0, noCoverages=0;
  unsigned int cov;

  bool wasNul=true;
  string::size_type prevNulpos=0;

  vector<unsigned int> covhisto;
  covhisto.resize(65535);

  for(string::size_type pos = 0; pos < d_mapping.size(); ++pos) {
    cov = d_mapping[pos].coverage;
    bool noCov = cov < 2;
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
        g_unm.push_back(unm);
      }
      wasNul=false;
    }
    else if(noCov && !wasNul) {
      wasNul=true;
      prevNulpos = pos;
    }
  }

  cerr << (boost::format("Average depth: %|40t|    %10.2f\n") % (1.0*totCoverage/d_mapping.size())).str();
  cerr << (boost::format("Undercovered nucleotides: %|40t| %10d (%.2f%%), %d ranges\n") % noCoverages % (noCoverages*100.0/d_mapping.size()) % g_unm.size()).str();

  uint64_t total = std::accumulate(covhisto.begin(), covhisto.end(), 0);

  fprintf(jsfp, "var histo=[");
  for(unsigned int i = 0; i < 250; i++ ) 
  {
    if(i)
      fputc(',', jsfp);
    fprintf(jsfp, "[%u, %f]", i, 1.0*covhisto[i]/total);
  }
  fprintf(jsfp,"];\n");
}


struct LociStats
{
  vector<boost::tuple<char,char,char>> samples;
};

typedef map<dnapos_t, LociStats> locimap_t;
locimap_t locimap;

// 0 if nothing interesting, positive if our read has insert at that position, negative if we have a delete at that position
int MBADiff(dnapos_t pos, const FastQRead& fqr, const string& reference)
{
  string::size_type n, m, d;
  int sn, i;
  struct varray *ses = varray_new(sizeof(struct diff_edit), NULL);
  
  n = reference.length();
  m = fqr.d_nucleotides.length();
  if ((d = diff(reference.c_str(), 0, n, fqr.d_nucleotides.c_str(), 0, m, NULL, NULL, NULL, 0, ses, &sn, NULL)) == -1) {
    MMNO(errno);
    printf("Error\n");
    return EXIT_FAILURE;
  }
  int ret=0;
  if(sn == 4 && d == 2) {
    struct diff_edit *match1 = (struct diff_edit*)varray_get(ses, 0),
      *change1=(struct diff_edit*)varray_get(ses, 1), 
      *match2=(struct diff_edit*)varray_get(ses, 2), 
      *change2=(struct diff_edit*)varray_get(ses, 3);
    
    if(match1->op == DIFF_MATCH && match2->op==DIFF_MATCH) {
      if(change1->op == DIFF_DELETE && change2->op==DIFF_INSERT) {
	//        cout << "Have delete of "<<change1->len<<" in our read at " << pos+change1->off <<endl;
        ret=-change1->off;
      }
      else if(change1->op == DIFF_INSERT && change2->op==DIFF_DELETE) {
	//        cout<<"Have insert of "<<change1->len<<" in our read at "<<pos+change1->off<<endl;
        ret=change1->off;
      }
    }
  }
 
  //  printf("pos %u, d=%lu sn=%d\nUS:  %s\nREF: %s\n", pos, d, 
  //     sn, fqr.d_nucleotides.c_str(), reference.c_str());

  if(sn > 6)
    return 0;


  // should have 4 components, a MATCH, an INSERT/DELETE followed by a MATCH, followed by a DELETE/INSERT (inverse)

  for (i = 0; i < sn; i++) {
    struct diff_edit *e = (struct diff_edit*)varray_get(ses, i);

    switch (e->op) {
    case DIFF_MATCH:
      //      printf("MAT: ");
      //      fwrite(fqr.d_nucleotides.c_str() + e->off, 1, e->len, stdout);
      break;
    case DIFF_INSERT:
      //      printf("INS: ");
      //      fwrite(reference.c_str() + e->off, 1, e->len, stdout);
      break;
    case DIFF_DELETE:
      //      printf("DEL: ");
      //      fwrite(fqr.d_nucleotides.c_str() + e->off, 1, e->len, stdout);
      break;
    }
    //    printf("\n");
  }
  varray_del(ses);
  return ret;
}

unsigned int diffScore(ReferenceGenome& rg, dnapos_t pos, FastQRead& fqfrag)
{
  unsigned int diffcount=0;
  string reference = rg.snippet(pos, pos + fqfrag.d_nucleotides.length());
  for(string::size_type i = 0; i < fqfrag.d_nucleotides.size() && i < reference.size();++i) {
    if(fqfrag.d_nucleotides[i] != reference[i] && fqfrag.d_quality[i]>'@') 
      diffcount++;
  }

  if(diffcount >= 5) { // bit too different, try mbadiff!
    int res=MBADiff(pos, fqfrag, reference);
    if(res < 0 || res > 0)
      return 1;
  }

  return diffcount;
}

vector<dnapos_t> getTriplets(const vector<pair<dnapos_t, char>>& together, unsigned int shift) 
{
  vector<dnapos_t> ret;
  bool doPrint=false;
  for(unsigned int i = 0; i < together.size() - 2; ++i) {
    if(together[i].second=='L' && together[i+1].second=='M' && together[i+2].second=='R' && 
       together[i+1].first - together[i].first < 60 && together[i+2].first - together[i+1].first < 60) {
      dnapos_t lpos;
      lpos=together[i].first; 

      if(lpos < shift)
        continue;
      
      ret.push_back(lpos-shift);
    }
  }         
  if(doPrint) {
    for(auto i : together) {
      cout<<i.first<<i.second<<" ";
    }
    cout<<endl;
    cout<<"Returning: ";
    for(auto i : ret) {
      cout<<i<<" ";
    }
    cout<<endl;
  }
  return ret;
}

map<dnapos_t, unsigned int> g_insertCounts;

void printCorrectMappings(FILE* jsfp, const ReferenceGenome& rg, const std::string& name)
{
  fprintf(jsfp, "var %s=[", name.c_str());
  for(unsigned int i=0; i< 150;++i) {
    double total=rg.d_correctMappings[i] + rg.d_wrongMappings[i];
    double error= rg.d_wrongMappings[i]/total;
    double qscore=-10*log(error)/log(10);
    fprintf(jsfp, "%s[%d,%.2f]", i ? "," : "", i, qscore);
    cout<<"total "<<total<<", error: "<<error<<", qscore: "<<qscore<<endl;
  }
  fprintf(jsfp,"];\n");
}

string DNADiff(ReferenceGenome& rg, dnapos_t pos, FastQRead& fqfrag)
{
  string reference = rg.snippet(pos, pos + fqfrag.d_nucleotides.length());

  int diffcount=0;
  for(string::size_type i = 0; i < fqfrag.d_nucleotides.size() && i < reference.size();++i) {
    if(fqfrag.d_nucleotides[i] != reference[i] && fqfrag.d_quality[i]>'@') 
      diffcount++;
  }

  if(diffcount < 5) 
    rg.mapFastQ(pos, fqfrag);
  else {
    int res=MBADiff(pos, fqfrag, reference);
    if(res) {
      rg.mapFastQ(pos, fqfrag, res);
      diffcount=1;
      if(res > 0) { // our read has an insert at this position
        fqfrag.d_nucleotides.erase(res, 1); // this makes things align again
        fqfrag.d_quality.erase(res, 1); 
        g_insertCounts[pos+res]++;
      } else {      // our read has an erase at this position
        fqfrag.d_nucleotides.insert(-res, 1, 'X');
        fqfrag.d_quality.insert(-res, 1, 'X');
      }
    }
  }

  string diff;
  diff.reserve(fqfrag.d_nucleotides.length());

  for(string::size_type i = 0; i < fqfrag.d_nucleotides.size() && i < reference.size();++i) {
    if(fqfrag.d_nucleotides[i] != reference[i]) {
      diff.append(1, fqfrag.d_quality[i] > '@' ? '!' : '^');
      if(fqfrag.d_quality[i]>'@' && diffcount < 5 && ( (!fqfrag.reversed && i > 14) || (fqfrag.reversed && i < (150-14) ) ) ) 
        locimap[pos+i].samples.push_back(boost::make_tuple(fqfrag.d_nucleotides[i], fqfrag.d_quality[i], fqfrag.reversed ^ (i>75))); // head or tail
      
      if(diffcount < 5)
	rg.d_wrongMappings[fqfrag.reversed ?  (149-i) : i]++;
    }
    else {
      diff.append(1, ' ');
      rg.cover(pos+i,fqfrag.d_quality[i]);
      if(diffcount < 5)
	rg.d_correctMappings[fqfrag.reversed ?  (149-i) : i]++;
    }
  }
  if(diffcount > 5) {
    cout<<"US:  "<<fqfrag.d_nucleotides<<endl<<"DIF: ";
    cout<<diff<<endl<<"REF: "<<reference<<endl;
    cout<<"QUA: "<<fqfrag.d_quality<<endl;
    cout<<endl;
  }
  return diff;
}


void printUnmatched(FILE*fp, ReferenceGenome& rg, FASTQReader& fastq, const string& name)
{
  unsigned int unmcount=0;
  for(auto& unm : g_unm) {
    fprintf(fp, "%s[%d]=[", name.c_str(), unmcount++);
    for(dnapos_t pos = unm.pos - 500; pos < unm.pos + 500; ++pos) {
      if(pos != unm.pos - 500) 
        fprintf(fp, ", ");
      fprintf(fp, "[%d, %d]", pos, rg.d_mapping[pos].coverage);
    }
    fprintf(fp,"];\n");
    //rg.printFastQs(unm.pos, fastq);
  }
  fputs("\n", fp);
  fflush(fp);
}

unsigned int variabilityCount(const ReferenceGenome& rg, dnapos_t position, const LociStats& lc, double* fraction)
{
  vector<int> counts(256);
  counts[rg.snippet(position, position+1)[0]]+=rg.d_mapping[position].coverage;
  
  int forwardCount=0;

  for(auto& j : lc.samples) {
    counts[j.get<0>()]++;
    if(j.get<2>())
      forwardCount++;
  }
  sort(counts.begin(), counts.end());
  unsigned int nonDom=0;
  for(unsigned int i=0; i < 255; ++i) {
    nonDom+=counts[i];
  }
  
  if(nonDom + counts[255] < 20)
    return 0;

  *fraction = 1.0*forwardCount / (1.0*lc.samples.size());
  if(*fraction < 0.05 || *fraction > 0.9)
    return 0;

  return 100*nonDom/counts[255];
}


dnapos_t halfFind(const std::vector<int64_t>& fqpositions, FASTQReader &fastq, ReferenceGenome& rg, int keylen)
{
  boost::progress_display fuzzyProgress(fqpositions.size(), cerr);
  vector<dnapos_t> stillUnfound;
  FastQRead fqfrag;
  for(auto fqpos : fqpositions) {
    ++fuzzyProgress;
    fastq.seek(fqpos);
    fastq.getRead(&fqfrag);

    if(!rg.getReadPositions(fqfrag.d_nucleotides.substr(0, keylen)).empty())
      continue;
    if(!rg.getReadPositions(fqfrag.d_nucleotides.substr(keylen, keylen)).empty())
      continue;

    fqfrag.reverse();
    if(!rg.getReadPositions(fqfrag.d_nucleotides.substr(0, keylen)).empty())
      continue;
    if(!rg.getReadPositions(fqfrag.d_nucleotides.substr(keylen, keylen)).empty())
      continue;
    
    stillUnfound.push_back(fqpos);
  }
  return fqpositions.size()  - stillUnfound.size();
}


int fuzzyFind(std::vector<uint64_t>* fqpositions, FASTQReader &fastq, ReferenceGenome& rg, int keylen)
{
  boost::progress_display fuzzyProgress(fqpositions->size(), cerr);
  vector<uint64_t> stillUnfound;
  FastQRead fqfrag;
  string left, middle, right;
  typedef pair<dnapos_t, char> tpos;
  vector<dnapos_t> lpositions, mpositions, rpositions;
  unsigned int fuzzyFound=0;
  for(auto fqpos : *fqpositions) {
    fastq.seek(fqpos);
    fastq.getRead(&fqfrag);
    
    struct Match 
    {
      Match() = default;
      Match(unsigned int score_, bool reversed_) : score(score_), reversed(reversed_){}
      unsigned int score;
      bool reversed;
    };
    map<dnapos_t, Match> allMatches;
    for(unsigned int attempts=0; attempts < 12; ++attempts) {
      for(int tries = 0; tries < 2; ++tries) {
        if(tries)
          fqfrag.reverse();
        left=fqfrag.d_nucleotides.substr(attempts*3, keylen);   middle=fqfrag.d_nucleotides.substr(50+attempts*3, keylen); right=fqfrag.d_nucleotides.substr(100+attempts*3, keylen);
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

        auto matches=getTriplets(together, 3*attempts);
        int score;
        for(auto match : matches) {
          if(allMatches.count(match))
            continue;
             
          score = diffScore(rg, match, fqfrag);
          
          allMatches[match] = Match{score, fqfrag.reversed};
          if(score==0) // won't get any better than this
            goto done;
        }
      }
    }
  done:;
    if(!allMatches.empty()) {
      //      cout<<"Have "<<allMatches.size()<<" different matches"<<endl;
      if(allMatches.size()==1) {
        // cout<<"  Position "<<allMatches.begin()->first<<" had score "<<allMatches.begin()->second.score<<endl;
        if(fqfrag.reversed != allMatches.begin()->second.reversed)
          fqfrag.reverse();

        DNADiff(rg, allMatches.begin()->first, fqfrag);
      } 
      else {
        map<unsigned int, vector<pair<dnapos_t, bool>>> scores;
        for(auto match: allMatches) {
	  //          cout<<"  Position "<<match.first<<" had score "<<match.second.score<<endl;
          scores[match.second.score].push_back(make_pair(match.first, match.second.reversed));
        }
        
        const auto& first = scores.begin()->second;
        auto pick = first[random() % first.size()];
	//        cout<<" Picking: "<< pick.first <<endl;
        if(fqfrag.reversed != pick.second)
          fqfrag.reverse();
        DNADiff(rg, pick.first, fqfrag);
      }
      fuzzyFound++;
    }
    else {
      stillUnfound.push_back(fqpos);
    }

    ++fuzzyProgress;
  }
  stillUnfound.swap(*fqpositions);
  cerr<<"\r";
  
  return fuzzyFound;
}

typedef vector<accumulator_set<double, stats<tag::mean, tag::moment<2> > > > qstats_t;

void printQualities(FILE* jsfp, const qstats_t& qstats)
{
  int i=0;
  fprintf(jsfp, "qualities=[");
  for(const auto& q : qstats) {
    if(i)
      fputs(",", jsfp);
    fprintf(jsfp, "[%d, %f]", i, mean(q));
    ++i;
  }
  fputs("];\n", jsfp);

  i=0;
  fprintf(jsfp, "qhilo=[");
  for(const auto& q : qstats) {
    if(i++)
      fputs(",", jsfp);
    fprintf(jsfp, "[%f, %f]", mean(q)-sqrt(moment<2>(q)), mean(q)+sqrt(moment<2>(q)));
  }
  fputs("];\n", jsfp);
  fflush(jsfp);
      
}

int main(int argc, char** argv)
{
  srandom(time(0));
  

  if(argc!=3) {
    fprintf(stderr, "Syntax: powerdna fastqfile fastafile\n");
    exit(EXIT_FAILURE);
  }
  GeneAnnotationReader gar("NC_012660.csv");
  cerr<<"Done reading "<<gar.size()<<" annotations"<<endl;
  
  FASTQReader fastq(argv[1]);

  unsigned int bytes=0;
  FastQRead fqfrag;
  bytes=fastq.getRead(&fqfrag); // get a read to index based on its size
  
  ReferenceGenome rg(argv[2]);

  rg.index(fqfrag.d_nucleotides.size());

  cerr<<"Loading positive control filter genome(s)"<<endl;
  ReferenceGenome phix("phi-x174.fasta");
  phix.index(fqfrag.d_nucleotides.size());

  dnapos_t pos;

  uint64_t withAny=0, found=0, notFound=0, total=0, qualityExcluded=0, fuzzyFound=0, phixFound=0;

  unique_ptr<FILE, int(*)(FILE*)> jsfp(fopen("data.js","w"), fclose);

  cerr<<"Performing exact matches of reads to reference genome"<<endl;
  boost::progress_display show_progress( filesize(argv[1]), cerr);
 
  qstats_t qstats;
  qstats.resize(fqfrag.d_nucleotides.size());
  accumulator_set<double, stats<tag::mean, tag::moment<2> > > qstat;
  vector<unsigned int> qcounts(256);
  vector<uint64_t> unfoundReads;
  do { 
    show_progress += bytes;
    total++;
    for(string::size_type pos = 0 ; pos < fqfrag.d_quality.size(); ++pos) {
      unsigned int i = fqfrag.d_quality[pos]-33;
      qstat(i);
      qstats[pos](i);
      qcounts[i]++;
    }
    
    if(fqfrag.d_nucleotides.find('N') != string::npos) {
      withAny++;
      continue;
    }
    pos = rg.getReadPosBoth(&fqfrag);
    if(pos == dnanpos ) {
      if(phix.getReadPosBoth(&fqfrag)!=dnanpos) {
        phixFound++;
      }
      else {
        unfoundReads.push_back(fqfrag.position);
        notFound++;
      }
    }
    else {
      rg.mapFastQ(pos, fqfrag);
      found++;
    }
  } while((bytes=fastq.getRead(&fqfrag)));

  uint64_t totNucleotides=total*fqfrag.d_nucleotides.length();
  fprintf(jsfp.get(), "qhisto=[");
  for(int c=0; c < 255; ++c) {
    fprintf(jsfp.get(), "%s[%d,%f]", c ? "," : "", (int)c, 1.0*qcounts[c]/totNucleotides);
  }
  fprintf(jsfp.get(),"];");


  cerr << (boost::format("Total reads: %|40t| %10d (%.2f gigabps)") % total % (totNucleotides/1000000000.0)).str() <<endl;
  cerr << (boost::format("Excluded control reads: %|40t|-%10d") % phixFound).str() <<endl;
  cerr << (boost::format("Quality excluded: %|40t|-%10d") % qualityExcluded).str() <<endl;
  cerr << (boost::format("Ignored reads with N: %|40t|-%10d") % withAny).str()<<endl;
  cerr << (boost::format("Full matches: %|40t|-%10d (%.02f%%)\n") % found % (100.0*found/total)).str();
  cerr << (boost::format("Not fully matched: %|40t|=%10d (%.02f%%)\n") % notFound % (notFound*100.0/total)).str();

  cerr << (boost::format("Mean Q: %|40t|    %10.2f +- %.2f\n") % mean(qstat) % sqrt(moment<2>(qstat))).str();

  for(auto& i : rg.d_correctMappings) {
    i=found;
  }

  for(auto& i : phix.d_correctMappings) {
    i=phixFound;
  }

  rg.printCoverage(jsfp.get());
  
  printQualities(jsfp.get(), qstats);
  printUnmatched(jsfp.get(), rg, fastq, "perfect");
  cout<<"---------"<<endl;

  int keylen=11;
  rg.index(keylen);
  phix.index(keylen);

  cerr<<"Performing sliding window partial matches"<<endl;

  fuzzyFound = fuzzyFind(&unfoundReads, fastq, rg, 11);
  auto fuzzyPhixFound = fuzzyFind(&unfoundReads, fastq, phix, 11);
  cerr<<fuzzyPhixFound<<" fuzzy @ phix"<<endl;

  cerr<<(boost::format("Fuzzy found: %|40t|-%10d\n")%fuzzyFound).str();  
  cerr<<(boost::format("Unmatchable reads:%|40t|=%10d (%.2f%%)\n") 
         % unfoundReads.size() % (100.0*unfoundReads.size()/total)).str();

  /*
  rg.index(75);
  auto chimericFound=halfFind(unfoundReads, fastq, rg, 75);
  cerr<<(boost::format("Probable chimeric reads:%|40t|-%10d (%.2f%%)\n") 
         % chimericFound % (100.0*chimericFound/total)).str();
  */

  /*
  rg.printFastQs(5718000, fastq);  
  */
  /*
  FILE *fp=fopen("unfound.fastq", "w");
  for(const auto& pos :  unfoundReads) {
    fastq.seek(pos);
    fastq.getRead(&fqfrag);
    fprintf(fp, "%s\n%s\n", fqfrag.d_nucleotides.c_str(), fqfrag.d_quality.c_str());
  }
  fclose(fp);
  */
  cout<<"After sliding matching: "<<endl;
  printUnmatched(jsfp.get(), rg, fastq, "fuzzy");

  printCorrectMappings(jsfp.get(), rg, "referenceQ");
  printCorrectMappings(jsfp.get(), phix, "controlQ");
  cerr<<"Found "<<locimap.size()<<" varying loci"<<endl;
  uint64_t seriouslyVariable=0;
  boost::format fmt1("%-10d: %3d*%c ");
  string fmt2("                  ");
  int aCount, cCount, tCount, gCount;
  double fraction;
  for(locimap_t::iterator iter = locimap.begin(); iter != locimap.end(); ++iter) {
    unsigned int varcount=variabilityCount(rg, iter->first, iter->second, &fraction);
    if(varcount < 20) 
      continue;
    char c=rg.snippet(iter->first, iter->first+1)[0];
    aCount = cCount = tCount = gCount = 0;
    switch(c) {
    case 'A':
      aCount += rg.d_mapping[iter->first].coverage;
      break;
    case 'C':
      cCount += rg.d_mapping[iter->first].coverage;
      break;
    case 'T':
      tCount += rg.d_mapping[iter->first].coverage;
      break;
    case 'G':
      gCount += rg.d_mapping[iter->first].coverage;
      break;
    }

    cout<< (fmt1 % iter->first % rg.d_mapping[iter->first].coverage % rg.snippet(iter->first, iter->first+1) ).str();
    sort(iter->second.samples.begin(), iter->second.samples.end());

    seriouslyVariable++;
    for(vector<boost::tuple<char,char,char> >::const_iterator j = iter->second.samples.begin(); 
        j != iter->second.samples.end(); ++j) {
      c=j->get<0>();
      switch(c) {
      case 'A':
        aCount++;
        break;
      case 'C':
        cCount++;
        break;
      case 'T':
        tCount++;
        break;
      case 'G':
        gCount++;
        break;
      }
      
      cout<<c;
    }
    cout<<endl<<fmt2;
    for(vector<boost::tuple<char,char,char> >::const_iterator j = iter->second.samples.begin(); 
        j != iter->second.samples.end(); ++j) {
      cout<<j->get<1>();
    }
    cout<<endl<<fmt2;
    for(vector<boost::tuple<char,char,char> >::const_iterator j = iter->second.samples.begin(); 
        j != iter->second.samples.end(); ++j) {
      cout<< (j->get<2>() ? 'R' : '.');
    }

    int tot=iter->second.samples.size() + rg.d_mapping[iter->first].coverage;
    cout<<endl;
    vector<GeneAnnotation> gas=gar.lookup(iter->first);
    if(!gas.empty()) {
      cout<<fmt2<<"Annotation: ";
      for(auto& ga : gas) {
        cout<<ga.name<<" ["<<ga.tag<<"], ";
      }
      cout<<endl;
    }
    cout<<fmt2<< "Fraction tail: "<<fraction<<", "<<iter->second.samples.size()<<endl;
    cout<<fmt2<< "A: " << aCount*100/tot <<"%, C: "<<cCount*100/tot<<"%, G: "<<gCount*100/tot<<"%, T: "<<tCount*100/tot<<"%"<<endl;

    rg.printFastQs(iter->first, fastq);
  }
  cerr<<"Found "<<seriouslyVariable<<" seriously variable loci"<<endl;
  cerr<<"Found "<<g_insertCounts.size()<<" loci with at least one insert in a read"<<endl;
  struct revsort
  {
    bool operator()(const unsigned int&a, const unsigned int&b) const
    { return a > b;} 
  };
  map<unsigned int, vector<dnapos_t>, revsort> topInserts;

  unsigned int seriousInserts=0;
  for(const auto& insloc : g_insertCounts) {
    topInserts[insloc.second].push_back(insloc.first);
    if(insloc.second > 10)
      seriousInserts++;
  }
  cerr<<"Found "<<seriousInserts<<" serious inserts"<<endl;
  
  for(const auto& insert : topInserts) {
    if(insert.first < 5)
      break;
    for(const auto& position : insert.second) {
      cout<<position<<"\t"<<insert.first<<" inserts"<<endl;
      rg.printFastQs(position, fastq);
    }
  }
  exit(EXIT_SUCCESS);
}
