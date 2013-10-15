#include <stdio.h>
#include <iostream>
#include <stdint.h>
#include <map>
#include <vector>
#include <string>
#include <string.h>
#include <stdexcept>
#include <boost/foreach.hpp>
#include <boost/progress.hpp>
#include <algorithm>
#include <numeric>
#include <unistd.h>
#include <math.h>
#include <boost/format.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/density.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include "geneannotated.hh"
#include "misc.hh"
#include "fastq.hh"

using namespace boost;
using namespace boost::accumulators;

extern "C" {
#include "hash.h"
}
using namespace std;

 
class ReferenceGenome
{
public:
  ReferenceGenome(FILE* fp);
  uint64_t size() const {
    return d_genome.size();
  }

  vector<uint64_t> getReadPositions(const std::string& nucleotides)
  {
    vector<uint64_t> ret;
    if(nucleotides.length() != d_indexlength)
      throw runtime_error("Attempting to find fragment of length we've not indexed for");
      
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
  
  uint64_t getReadPosBoth(const FastQRead& fqp) // tries both
  {
    FastQRead fq = fqp;

    vector<uint64_t> positions;
    string nucleotides;
    for(int tries = 0; tries < 2; ++tries) {
      positions = getReadPositions(fq.d_nucleotides);

      if(!positions.empty()) {
	if(positions.size() > 1)
	  random_shuffle(positions.begin(), positions.end());

	cover(positions[0], fq.d_nucleotides.size(), fq.d_quality);
	return positions[0];
      }
      fq.reverse();
    }
    return string::npos;
  }

  void cover(uint64_t pos, unsigned int length, const std::string& quality) 
  {
    const char* p = quality.c_str();
    for(unsigned int i = 0; i < length; ++i) {
      if(p[i]-33 > 30)
	d_coverage[pos+i]++;
    }
  }

  void cover(uint64_t pos, char quality) 
  {
    if(quality-33 > 30)
      d_coverage[pos]++;
  }


  string snippet(uint64_t start, uint64_t stop) const { 
    return d_genome.substr(start, stop-start);
  }
  void printCoverage();
  void index(int length);

  vector<uint16_t> d_coverage;
private:
  unsigned int d_indexlength;
  string d_genome;

  struct HashPos {
    HashPos(uint32_t hash_, uint64_t pos) : d_hash(hash_), d_pos(pos)
    {}
    HashPos(){}
    uint32_t d_hash;
    uint64_t d_pos;
    
    bool operator<(const HashPos& rhs) const 
    {
      return d_hash < rhs.d_hash;
    }
  };

  typedef vector<HashPos> index_t;
  index_t d_index;
};

ReferenceGenome::ReferenceGenome(FILE *fp)
{
  char line[256]="";

  sfgets(line, sizeof(line), fp);
  if(line[0] != '>') 
    throw runtime_error("Input not FASTA");

  while(fgets(line, sizeof(line), fp)) {
    chomp(line);
    d_genome.append(line);
  }
  d_coverage.resize(d_genome.size());
  d_indexlength=0;
}

void ReferenceGenome::index(int length)
{
  d_index.clear();
  d_indexlength=length;
  cerr<<"Indexing "<<d_genome.length()<<" nucleotides for length "<<length<<endl;
  boost::progress_display show_progress( d_genome.length() - d_indexlength, cerr);

  for(string::size_type pos = 0 ; pos < d_genome.length() - length; ++pos) {
    uint32_t hashval = hash(d_genome.c_str() + pos, d_indexlength, 0);
    d_index.push_back(HashPos(hashval, pos));
    ++show_progress;
  }

  cerr<<endl<<"Sorting hashes..";
  sort(d_index.begin(), d_index.end());
  cerr<<" done"<<endl;
    
  cerr<<"Average hash fill: "<<1.0*d_genome.length()/d_index.size()<<endl;
}



struct Rematched
{
  string left, middle, right;
};

struct Unmatched
{
  string left, unmatched, right;
  uint64_t pos;
  vector<Rematched> matches;
};

vector<Unmatched> g_unm;

void ReferenceGenome::printCoverage()
{
  uint64_t totCoverage=0, noCoverages=0;
  unsigned int cov;

  vector<unsigned int> nulls;
  nulls.resize(2000);
  unsigned int binwidth = d_coverage.size()/2000;

  bool wasNul=true;
  string::size_type prevNulpos=0;

  vector<unsigned int> covhisto;
  covhisto.resize(65535);

  for(string::size_type pos = 0; pos < d_coverage.size(); ++pos) {
    cov = d_coverage[pos];
    bool noCov = cov < 1;
    covhisto[cov]++;
    totCoverage += cov;

    if(noCov) {
      noCoverages++;
      nulls[pos / binwidth]++;
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

  cerr<<"Average depth: "<<totCoverage/d_coverage.size()<<endl;
  cerr<<"No coverage: "<<noCoverages*100.0/d_coverage.size()<<"%"<<endl;

  uint64_t total = std::accumulate(covhisto.begin(), covhisto.end(), 0);

  cout<<"var histo=[";
  for(unsigned int i = 0; i < 250; i++ ) 
  {
    if(i)
      cout<<',';
    std::cout << '['<< i << ',' << 1.0*covhisto[i]/total << ']';
  }
  cout<<"];\n";

  //  for(unsigned int n = 0; n < 2000; ++n) 
  //  cout<<n*binwidth<<"\t"<<nulls[n]<<endl;  
}


struct LociStats
{
  vector<pair<char,char > > samples;
};

typedef map<uint64_t, LociStats> locimap_t;
locimap_t locimap;

string DNADiff(ReferenceGenome& rg, uint64_t pos, const std::string& a, const std::string& b, const std::string& quality)
{
  string diff;
  diff.reserve(a.length());

  int diffcount=0;
  for(string::size_type i = 0; i < a.size() && i < b.size();++i) {
    if(a[i] != b[i] && quality[i]>'@') 
      diffcount++;
  }
  cout<<"US:  "<<a<<endl<<"DIF: ";
  for(string::size_type i = 0; i < a.size() && i < b.size();++i) {
    if(a[i] != b[i]) {
      diff.append(1, quality[i] > '@' ? '!' : '^');
      if(quality[i]>'@' && diffcount < 5)
	locimap[pos+i].samples.push_back(make_pair(a[i], quality[i]));
    }
    else {
      diff.append(1, ' ');
      rg.cover(pos+i,quality[i]);
    }
  }

  cout<<diff<<endl<<"REF: "<<b<<endl;
  cout<<"QUA: "<<quality<<endl;
  return diff;
}


void printUnmatched(ReferenceGenome& rg, const string& name)
{
  unsigned int unmcount=0;
  BOOST_FOREACH(Unmatched& unm, g_unm) {
    printf("%s[%d]=[", name.c_str(), unmcount++);
    for(uint64_t pos = unm.pos - 500; pos < unm.pos + 500; ++pos) {
      if(pos != unm.pos - 500) 
	printf(", ");
      printf("[%ld, %d]", pos, rg.d_coverage[pos]);
    }
    printf("];\n");
  }
  cout<<endl;
}

unsigned int variabilityCount(const ReferenceGenome& rg, uint64_t position, const LociStats& lc)
{
  vector<int> counts(256);
  counts[rg.snippet(position, position+1)[0]]+=rg.d_coverage[position];

  for(vector<pair<char,char> >::const_iterator j = lc.samples.begin(); j!= lc.samples.end(); ++j) {
    counts[j->first]++;
  }
  sort(counts.begin(), counts.end());
  unsigned int nonDom=0;
  for(unsigned int i=0; i < 255; ++i) {
    nonDom+=counts[i];
  }
  
  return 100*nonDom/counts[255];

}

int main(int argc, char** argv)
{
  if(argc!=3) {
    fprintf(stderr, "Syntax: powerdna fastqfile fastafile\n");
    exit(EXIT_FAILURE);
  }
  GeneAnnotationReader gar("NC_012660.csv");
  cerr<<"Done reading "<<gar.size()<<" annotations"<<endl;
  
  FASTQReader fastq(argv[1]);
  FILE* fasta = fopen(argv[2], "r");

  ReferenceGenome rg(fasta);

  FILE* phixFP = fopen("phi-x174.fasta", "r");
  ReferenceGenome phix(phixFP);

  FastQRead fqfrag;
  
  unsigned int bytes=0;
  bytes=fastq.getRead(&fqfrag); // get a fragment to index based on its size
  rg.index(fqfrag.d_nucleotides.size());
  phix.index(fqfrag.d_nucleotides.size());

  uint64_t pos;

  uint64_t withAny=0, found=0, notFound=0, total=0, qualityExcluded=0, fuzzyFound=0, phixFound=0;
  boost::progress_display show_progress( filesize(argv[1]), cerr);

  accumulator_set<double, stats<tag::mean, tag::median > > acc;

  vector<uint64_t> unfoundReads;
  do { 
    show_progress += bytes;
    total++;
    BOOST_FOREACH(char c, fqfrag.d_quality) {
      double i = c-33;
      acc(i);
    }

    if(fqfrag.d_nucleotides.find('N') != string::npos) {
      withAny++;
      continue;
    }
    pos = rg.getReadPosBoth(fqfrag);
    if(pos == string::npos ) {
      if(phix.getReadPosBoth(fqfrag)!=string::npos) {
	phixFound++;
      }
      else {
	unfoundReads.push_back(fastq.getPos());
	notFound++;
      }
    }
    else {
      found++;
    }
  } while((bytes=fastq.getRead(&fqfrag)));
  cerr<<"Total fragments: "<<total<<endl;
  cerr<<"Phix: "<<phixFound<<endl;
  cerr<<"Quality Excluded: "<<qualityExcluded<<endl;
  cerr<<"Full matches: "<< found<< endl;
  cerr<< (boost::format("Not found: %d (%.02f%%)\n") % notFound % (notFound*100.0/total)).str() << endl;
  cerr<<"Reads with N: "<<withAny++<<endl;
  cerr<<"Not found: "<<unfoundReads.size()<<endl;
  cerr << "Mean Q:   " << mean(acc) << std::endl;
  cerr << "Median Q: " << median(acc) << std::endl;
  
  rg.printCoverage();
  
  printUnmatched(rg,"perfect");
  cout<<"---------"<<endl;

  int keylen=11;
  rg.index(keylen);

  string left, middle, right;
  vector<uint64_t> lpositions, mpositions, rpositions;
  vector<uint64_t> stillUnfound;
  typedef pair<uint64_t, char> tpos;

  boost::progress_display fuzzyProgress(unfoundReads.size(), cerr);

  BOOST_FOREACH(uint64_t pos, unfoundReads) {
    fastq.seek(pos);
    fastq.getRead(&fqfrag);

    for(unsigned int attempts=0; attempts < 12; ++attempts) {
      for(int tries = 0; tries < 2; ++tries) {
	if(tries)
	  fqfrag.reverse();
	left=fqfrag.d_nucleotides.substr(attempts*3, keylen);	middle=fqfrag.d_nucleotides.substr(50+attempts*3, keylen); right=fqfrag.d_nucleotides.substr(100+attempts*3, keylen);
      	lpositions=rg.getReadPositions(left); mpositions=rg.getReadPositions(middle); rpositions=rg.getReadPositions(right);

       	vector<tpos> together;
	BOOST_FOREACH(uint64_t fpos, lpositions) {   together.push_back(make_pair(fpos, 'L')); }	
	BOOST_FOREACH(uint64_t fpos, mpositions) {   together.push_back(make_pair(fpos, 'M')); }
	BOOST_FOREACH(uint64_t fpos, rpositions) {   together.push_back(make_pair(fpos, 'R')); }
	if(together.size() < 3) 
	  continue;
	
	sort(together.begin(), together.end());

	for(unsigned int i = 0; i < together.size() - 2; ++i) {
	  if(together[i].second=='L' && together[i+1].second=='M' && together[i+2].second=='R' && 
	     together[i+1].first - together[i].first < 60 && together[i+2].first - together[i+1].first < 60) {
	    uint64_t lpos, mpos, rpos;
	    lpos=together[i].first; 	    mpos=together[i+1].first; 	    rpos=together[i+2].first;
	    cout<<lpos<<together[i].second<<" - " << mpos<<together[i+1].second<<
	      " - " << rpos<<together[i+2].second<< " (" << 
	      mpos - lpos<< " - " << rpos - mpos <<"), reversed: "<<tries<<endl;
	    
	    if(lpos < 3*attempts)
	      continue;
	    
	    uint64_t matchOffset=lpos - 3*attempts;
	    DNADiff(rg, matchOffset, fqfrag.d_nucleotides, rg.snippet(matchOffset, matchOffset+150), fqfrag.d_quality);
	    fuzzyFound++;
	    goto foundIt;
	    break;
	  }
	}
      }
    }
    stillUnfound.push_back(pos);
  foundIt:;
    ++fuzzyProgress;
  }
  cerr<<"\rFuzzy found: "<<fuzzyFound<<endl;
  fuzzyFound=0;
  unfoundReads.swap(stillUnfound);
  cerr<<"Have "<<unfoundReads.size()<<" unfound reads left"<<endl;
  
  FILE *fp=fopen("unfound.fastq", "w");
  BOOST_FOREACH(uint64_t pos, unfoundReads) {
    fastq.seek(pos);
    fastq.getRead(&fqfrag);
    fprintf(fp, "%s\n%s\n", fqfrag.d_nucleotides.c_str(), fqfrag.d_quality.c_str());
  }
  fclose(fp);
  cout<<"After sliding matching: "<<endl;
  printUnmatched(rg, "fuzzy");

  cerr<<"Found "<<locimap.size()<<" varying loci"<<endl;
  uint64_t seriouslyVariable=0;
  boost::format fmt1("%-10d: %3d*%c ");
  string fmt2("                  ");
  int aCount, cCount, tCount, gCount;

  for(locimap_t::iterator iter = locimap.begin(); iter != locimap.end(); ++iter) {
    unsigned int varcount=variabilityCount(rg, iter->first, iter->second);
    if(varcount < 20) 
      continue;
    char c=rg.snippet(iter->first, iter->first+1)[0];
    aCount = cCount = tCount = gCount = 0;
    switch(c) {
    case 'A':
      aCount += rg.d_coverage[iter->first];
      break;
    case 'C':
      cCount += rg.d_coverage[iter->first];
      break;
    case 'T':
      tCount += rg.d_coverage[iter->first];
      break;
    case 'G':
      gCount += rg.d_coverage[iter->first];
      break;
    }

    cout<< (fmt1 % iter->first % rg.d_coverage[iter->first] % rg.snippet(iter->first, iter->first+1) ).str();
    sort(iter->second.samples.begin(), iter->second.samples.end());

    seriouslyVariable++;
    for(vector<pair<char,char> >::const_iterator j = iter->second.samples.begin(); 
	j != iter->second.samples.end(); ++j) {
      c=j->first;
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
      
      cout<<j->first;
    }
    cout<<endl<<fmt2;
    for(vector<pair<char,char> >::const_iterator j = iter->second.samples.begin(); 
	j != iter->second.samples.end(); ++j) {
      cout<<j->second;
    }
    int tot=iter->second.samples.size() + rg.d_coverage[iter->first];
    cout<<endl;
    vector<GeneAnnotation> gas=gar.lookup(iter->first);
    if(!gas.empty()) {
      cout<<fmt2<<"Annotation: ";
      BOOST_FOREACH(const GeneAnnotation& ga, gas) {
	cout<<ga.name<<" ["<<ga.tag<<"], ";
      }
      cout<<endl;
    }
    
    cout<<fmt2<< "A: " << aCount*100/tot <<"%, C: "<<cCount*100/tot<<"%, G: "<<gCount*100/tot<<"%, T: "<<tCount*100/tot<<"%";

#if 0
    vector<int> suf;
    suf.push_back(100*aCount/tot);     suf.push_back(100*cCount/tot);    suf.push_back(100*gCount/tot);    
    suf.push_back(100*tCount/tot);
    sort(suf.begin(), suf.end());
    cout<<endl<<fmt2<<"PRO: ";
    BOOST_FOREACH(int k, suf) { cout << k<<" ";}

#endif
    cout<<endl;
  }
  cerr<<"Found "<<seriouslyVariable<<" seriously variable loci"<<endl;
  exit(EXIT_SUCCESS);

  //fastq.seek(0);  
#if 0
  string::size_type lpos, rpos;
  while((bytes=fastq.getRead(&fqfrag))) {
    BOOST_FOREACH(Unmatched& unm, g_unm) {
      lpos = rpos = string::npos;
      lpos = fqfrag.d_nucleotides.find(unm.left);
      if(lpos != string::npos)
	rpos = fqfrag.d_nucleotides.find(unm.right);
      
      if(rpos == string::npos) {
	fqfrag.reverse();
	lpos = fqfrag.d_nucleotides.find(unm.left);
	if(lpos != string::npos)
	  rpos = fqfrag.d_nucleotides.find(unm.right);
      }
      
      if(rpos == string::npos)
	continue;

      // lpos & rpos now correct 
      cout<<"UNM: "<<unm.pos << " - " << unm.pos + unm.unmatched.length() << " (len = "<< unm.unmatched.length()<<", lpos = "<<lpos<<", rpos = "<<rpos<<")"<<endl;
      cout << unm.left << " | "<< unm.unmatched<< " | " << unm.right<<endl;

      Rematched rem;
      if(lpos + 40 > rpos) {
	rem.left = fqfrag.d_nucleotides.substr(lpos, 40);  
	rem.middle="1";
	rem.right=fqfrag.d_nucleotides.substr(rpos);
      }
      else {
	rem.left = fqfrag.d_nucleotides.substr(lpos, 40);
	rem.middle = fqfrag.d_nucleotides.substr(lpos+40, rpos-lpos-40);
	rem.right = fqfrag.d_nucleotides.substr(rpos);
      }
      cout << rem.left << " | " << rem.middle << " | " << rem.right <<endl;
      unm.matches.push_back(rem);
    }
  }

  cout<<"------------------"<<endl;
  BOOST_FOREACH(Unmatched& unm, g_unm) {
    cout<<"REM: "<<unm.pos << " - " << unm.pos + unm.unmatched.length() << " (len = "<< unm.unmatched.length()<<", matches: "<<unm.matches.size()<<")\n";
    cout << "ref "<<unm.left << " | "<< unm.unmatched<< " | " << unm.right<<endl;

    BOOST_FOREACH(Rematched& rem, unm.matches) {
      cout <<"us  "<< rem.left << " | "<< rem.middle<< " | " << rem.right<<endl;
    }
    cout<<endl;
    
  }
#endif
}
