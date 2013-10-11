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
#include <sys/types.h>
#include <sys/stat.h>
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
using namespace boost;
using namespace boost::accumulators;


extern "C" {
#include "hash.h"
}
using namespace std;

struct FastQFragment
{
  string d_nucleotides;
  string d_quality;
  bool exceedsQuality(unsigned int);
  void reverse();
};

bool FastQFragment::exceedsQuality(unsigned int limit)
{
  uint8_t q;
  for(string::size_type pos = 0 ; pos < d_quality.size(); ++pos) {
    q = d_quality[pos] - 33;
    if(q < limit)
      return false;
  }
  return true;
}

void FastQFragment::reverse()
{
  std::reverse(d_nucleotides.begin(), d_nucleotides.end());
  for(string::iterator iter = d_nucleotides.begin(); iter != d_nucleotides.end(); ++iter) {
    if(*iter == 'C')
      *iter = 'G';
    else if(*iter == 'G')
      *iter = 'C';
    else if(*iter == 'A')
      *iter = 'T';
    else if(*iter == 'T')
      *iter = 'A';
  }
}

char* sfgets(char* p, int num, FILE* fp)
{
  char *ret = fgets(p, num, fp);
  if(!ret) 
    throw runtime_error("Unexpected EOF");
  return ret;
}
 
class ReferenceGenome
{
public:
  ReferenceGenome(FILE* fp);
  uint64_t size() const {
    return d_genome.size();
  }

  vector<uint64_t> getFragmentPositions(const std::string& nucleotides)
  {
    vector<uint64_t> ret;
    if(nucleotides.length() != d_indexlength)
      throw runtime_error("Attempting to find fragment of length we've not indexed for");
      
    uint32_t hashval = hash(nucleotides.c_str(), nucleotides.length(), 0);
    index_t::const_iterator iter = d_index.find(hashval);
    if(iter == d_index.end())
      return ret;

    BOOST_FOREACH(uint64_t off, iter->second) {
      if(!memcmp(d_genome.c_str() + off, nucleotides.c_str(), nucleotides.length())) {
	ret.push_back(off);
      }
    }
    return ret;
  }
  
  uint64_t getFragmentPos(const FastQFragment& fq, bool doCover=1)
  {
    vector<uint64_t> positions = getFragmentPositions(fq.d_nucleotides);
    if(doCover) {
      for(vector<uint64_t>::const_iterator iter = positions.begin(); iter != positions.end(); ++iter) {
	cover(*iter, fq.d_nucleotides.size(), fq.d_quality);
      }
    }
    if(positions.empty())
      return string::npos;
    return positions[0];
  }

  void cover(uint64_t pos, unsigned int length, const std::string& quality) 
  {
    const char* p = quality.c_str();
    for(unsigned int i = 0; i < length; ++i) {
      if(p[i]-33 > 30)
	d_coverage[pos+i]++;
    }
  }
  string snippet(uint64_t start, uint64_t stop) {
    return d_genome.substr(start, stop-start);
  }
  void printCoverage();
  void index(int length);

  vector<uint16_t> d_coverage;
private:
  unsigned int d_indexlength;
  string d_genome;

  typedef map<uint32_t, vector<uint64_t> > index_t;
  index_t d_index;
};

void chomp(char* line)
{
  char *p;
  p = strchr(line, '\r');
  if(p)*p=0;
  p = strchr(line, '\n');
  if(p)*p=0;
}

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
    d_index[hashval].push_back(pos);
    ++show_progress;
  }
  
  cerr<<endl;
  cerr<<"Average hash fill: "<<1.0*d_genome.length()/d_index.size()<<endl;
}

class FASTQReader
{
public:
  FASTQReader(const std::string& str) 
  {
    d_fp = fopen(str.c_str(), "r");
    if(!d_fp) 
      throw runtime_error("Unable to open file "+str+" for FASTQ inpuot");
  }

  uint64_t getPos() const
  {
    return d_pos;
  }

  void seek(uint64_t pos) 
  {
    fseek(d_fp, pos, SEEK_SET);
    // d_pos gets reset AFTER a read
  }

  unsigned int getFragment(FastQFragment* fq, unsigned int size=0);
private:
  FILE *d_fp;
  uint64_t d_pos;
};

unsigned int FASTQReader::getFragment(FastQFragment* fq, unsigned int size)
{
  d_pos = ftell(d_fp);
  char line[256]="";
  if(!fgets(line, sizeof(line), d_fp)) 
    return 0;
  if(line[0] != '@')
    throw runtime_error("Input not FASTQ");

  sfgets(line, sizeof(line), d_fp);
  chomp(line);
  
  if(size)
    fq->d_nucleotides.assign(line, line+size);
  else
    fq->d_nucleotides.assign(line);
  sfgets(line, sizeof(line), d_fp);
  sfgets(line, sizeof(line), d_fp);
  chomp(line);
  fq->d_quality.assign(line);
  return ftell(d_fp) - d_pos;
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

  cout<<"var histo=["<<endl;
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

uint64_t filesize(const char* name)
{
  struct stat buf;
  if(!stat(name, &buf)) {
    return buf.st_size;
  }
  return 0;
}


struct LociStats
{
  vector<pair<char,char > > samples;
};

typedef map<uint64_t, LociStats> locimap_t;
locimap_t locimap;

string DNADiff(uint64_t pos, const std::string& a, const std::string& b, const std::string& quality)
{
  string diff;
  diff.reserve(a.length());

  int diffcount=0;
  for(string::size_type i = 0; i < a.size() && i < b.size();++i) {
    if(a[i] != b[i]) 
      diffcount++;
  }

  for(string::size_type i = 0; i < a.size() && i < b.size();++i) {
    if(a[i] != b[i]) {
      diff.append(1, quality[i] > '@' ? '!' : '^');
      if(quality[i]>'@' && diffcount < 5)
	locimap[pos+i].samples.push_back(make_pair(a[i], quality[i]));
    }
    else
      diff.append(1, ' ');
  }
  return diff;
}


int main(int argc, char** argv)
{
  if(argc!=3) {
    fprintf(stderr, "Syntax: powerdna fastqfile fastafile\n");
    exit(EXIT_FAILURE);
  }
  GeneAnnotationReader gar("NC_012660.csv");
  cerr<<"Done reading annotations"<<endl;
  gar.lookup(10);
  vector<GeneAnnotation> ga = gar.lookup(1488618);
  cout<<ga[0].name<<endl;
  
  FASTQReader fastq(argv[1]);
  FILE* fasta = fopen(argv[2], "r");

  ReferenceGenome rg(fasta);
  FastQFragment fqfrag;
  
  unsigned int bytes=0;
  bytes=fastq.getFragment(&fqfrag); // get a fragment to index based on its size
  rg.index(fqfrag.d_nucleotides.size());

  uint64_t pos;
  uint64_t withAny=0, found=0, notFound=0, total=0, reverseFound=0, qualityExcluded=0, fuzzyFound=0;
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
    pos = rg.getFragmentPos(fqfrag);
    if(pos == string::npos) {
      fqfrag.reverse();
      pos = rg.getFragmentPos(fqfrag);
      if(pos != string::npos)
	reverseFound++;
      else {
	unfoundReads.push_back(fastq.getPos());
	notFound++;
      }
    }
    else {
      found++;
    }
  } while((bytes=fastq.getFragment(&fqfrag)));
  cerr<<"Total fragments: "<<total<<endl;
  cerr<<"Quality Excluded: "<<qualityExcluded<<endl;
  cerr<<"Found + Reverse found : "<< found<< " + " << reverseFound<< " = " <<found+reverseFound<<endl;
  cerr<< (boost::format("Not found: %d (%.02f%%)\n") % notFound % (notFound*100.0/total)).str() << endl;
  cerr<<"Fragments with N: "<<withAny++<<endl;
  std::cerr << "Mean:   " << mean(acc) << std::endl;
  std::cerr << "Median: " << median(acc) << std::endl;

  rg.printCoverage();
  
  unsigned int unmcount=0;
  BOOST_FOREACH(Unmatched& unm, g_unm) {
    printf("unm[%d]=[", unmcount++);
    for(uint64_t pos = unm.pos - 500; pos < unm.pos + 500; ++pos) {
      if(pos != unm.pos - 500) 
	printf(", ");
      printf("[%ld, %d]", pos, rg.d_coverage[pos]);
    }
    printf("];\n");
  }
  cout<<endl;
  cout<<"---------"<<endl;

  int keylen=11;
  rg.index(keylen);

  string left, middle, right;
  vector<uint64_t> lpos, mpos, rpos;
  BOOST_FOREACH(uint64_t pos, unfoundReads) {
    fastq.seek(pos);
    fastq.getFragment(&fqfrag);

    for(int tries = 0; tries < 2; ++tries) {
      left=fqfrag.d_nucleotides.substr(0, keylen);
      middle=fqfrag.d_nucleotides.substr(50, keylen);
      right=fqfrag.d_nucleotides.substr(100, keylen);

      lpos=rg.getFragmentPositions(left); 
      mpos=rg.getFragmentPositions(middle);
      rpos=rg.getFragmentPositions(right);

      //    cout << lpos.size()<<", "<<mpos.size()<<", " <<rpos.size()<<endl;
      typedef pair<uint64_t, char> tpos;

      vector<tpos> together;
      BOOST_FOREACH(uint64_t fpos, lpos) { 
	together.push_back(make_pair(fpos, 'L'));
      }
      BOOST_FOREACH(uint64_t fpos, mpos) { 
	together.push_back(make_pair(fpos, 'M'));
      }
      BOOST_FOREACH(uint64_t fpos, rpos) { 
	together.push_back(make_pair(fpos, 'R'));
      }
      sort(together.begin(), together.end());
      bool found=false;
      if(together.size()>=3) {
	for(unsigned int i = 0; i < together.size() - 3; ++i) {
	  if(together[i].second=='L' && together[i+1].second=='M' && together[i+2].second=='R' && together[i+1].first - together[i].first < 200 && together[i+2].first - together[i+1].first < 200) {
	    cout<<together[i].first<<together[i].second<<" - " <<
	      together[i+1].first<<together[i+1].second<<" - " <<
	      together[i+2].first<<together[i+2].second<<endl;
	    cout<<"REF: "<<rg.snippet(together[i].first, together[i+2].first+30)<<endl;
	    cout<<"US:  "<<fqfrag.d_nucleotides<<endl;
	    cout<<"DIF: "<<DNADiff(together[i].first, fqfrag.d_nucleotides, rg.snippet(together[i].first, together[i+2].first+30), fqfrag.d_quality)<<endl;
	    cout<<"QUA: "<<fqfrag.d_quality<<endl;
	    found=true;
	    fuzzyFound++;
	  }
	}
      }
      if(found)
	break;
      if(!tries)
	fqfrag.reverse();
      
      cout<<"UNF: ";
      BOOST_FOREACH(tpos& tp, together) {
	cout<<tp.first<<tp.second<<" ";
      }
      cout<<endl;
      
    }
  }
  cerr<<"Fuzzy found: "<<fuzzyFound<<endl;

  cerr<<"Found "<<locimap.size()<<" varying loci"<<endl;
  uint64_t seriouslyVariable=0;
  boost::format fmt1("%-10d: %3d*%c ");
  string fmt2("                  ");
  int aCount, cCount, tCount, gCount;

  
  for(locimap_t::iterator iter = locimap.begin(); iter != locimap.end(); ++iter) {
    if(rg.d_coverage[iter->first] < 10 ||  iter->second.samples.size() < 7) 
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
  while((bytes=fastq.getFragment(&fqfrag))) {
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
