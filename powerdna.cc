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

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/density.hpp>
#include <boost/accumulators/statistics/stats.hpp>

using namespace boost;
using namespace boost::accumulators;

typedef accumulator_set<unsigned int, features<tag::density> > acc;
typedef iterator_range<std::vector<std::pair<double, double> >::iterator > histogram_type; 
  


extern "C" {
#include "hash.h"
}
using namespace std;

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
  
  uint64_t getFragmentPos(const std::string& str) 
  {
    if(str.length() != d_indexlength)
      throw runtime_error("Attempting to find fragment of length we've not indexed for");
      
    uint32_t hashval = hash(str.c_str(), str.length(), 0);
    index_t::const_iterator iter = d_index.find(hashval);
    if(iter == d_index.end())
      return string::npos;

    const char* pos=0;
    BOOST_FOREACH(uint64_t off, iter->second) {
      if(!memcmp(d_genome.c_str() + off, str.c_str(), str.length())) {
	pos = d_genome.c_str() + off;
	cover(off, str.size());
      }
    }
    if(pos)
      return pos - d_genome.c_str();

    return string::npos;
  }

  void cover(uint64_t pos, unsigned int length) {
    while(length-- > 0) 
      d_coverage[pos++]++;
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

struct FastQFragment
{
  string d_nucleotides;
  void reverse();
};

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

unsigned int getFragment(FILE* fastq, FastQFragment* fq, unsigned int size=0)
{
  uint64_t pos = ftell(fastq);
  char line[256]="";
  if(!fgets(line, sizeof(line), fastq)) 
    return 0;
  if(line[0] != '@')
    throw runtime_error("Input not FASTQ");

  sfgets(line, sizeof(line), fastq);
  chomp(line);
  
  if(size)
    fq->d_nucleotides.assign(line, line+size);
  else
    fq->d_nucleotides.assign(line);
  sfgets(line, sizeof(line), fastq);
  sfgets(line, sizeof(line), fastq);
  return ftell(fastq) - pos;
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
  uint64_t totCoverage=0, noCoverage=0;
  unsigned int cov;

  vector<unsigned int> nulls;
  nulls.resize(2000);
  unsigned int binwidth = d_coverage.size()/2000;

  bool wasNul=true;
  string::size_type prevNulpos=0;

  vector<unsigned int> covhisto;
  covhisto.resize(65535);
  acc dens( tag::density::num_bins = 150, tag::density::cache_size = 10000);

  for(string::size_type pos = 0; pos < d_coverage.size(); ++pos) {
    cov = d_coverage[pos];
    covhisto[cov]++;
    totCoverage += cov;

    if(cov < 2) {
      noCoverage++;
      nulls[pos / binwidth]++;
    }
    
    if(cov && wasNul) {

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
    else if(!cov && !wasNul) {
      wasNul=true;
      prevNulpos = pos;
    }

  }

  double mu = totCoverage/d_coverage.size();
  cerr<<"Average depth: "<<totCoverage/d_coverage.size()<<endl;
  cerr<<"No coverage: "<<noCoverage*100.0/d_coverage.size()<<"%"<<endl;


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

int main(int argc, char** argv)
{
  if(argc!=3) {
    fprintf(stderr, "Syntax: powerdna fastqfile fastafile\n");
    exit(EXIT_FAILURE);
  }
  FILE* fastq = fopen(argv[1], "r");
  FILE* fasta = fopen(argv[2], "r");

  ReferenceGenome rg(fasta);
  FastQFragment fq;
  
  unsigned int bytes=0;
  bytes=getFragment(fastq, &fq); // get a fragment and index based on its size
  rg.index(fq.d_nucleotides.size());

  uint64_t pos;
  uint64_t withAny=0, found=0, notFound=0, total=0, reverseFound=0;
  boost::progress_display show_progress( filesize(argv[1]), cerr);

  do { 
    show_progress += bytes;
    total++;
    if(fq.d_nucleotides.find('N') != string::npos) {
      withAny++;
      continue;
    }
    pos = rg.getFragmentPos(fq.d_nucleotides);
    if(pos == string::npos) {
      fq.reverse();
      pos = rg.getFragmentPos(fq.d_nucleotides);
      if(pos != string::npos)
	reverseFound++;
      else
	notFound++;
    }
    else {
      found++;
    }
  } while((bytes=getFragment(fastq, &fq)));
  cerr<<"Total fragments: "<<total<<endl;
  cerr<<"Found: "<<found<<endl;
  cerr<<"Reverse found: "<<reverseFound<<endl;
  cerr<<"Not found: "<< notFound<<endl;
  cerr<<"Fragments with N: "<<withAny++<<endl;

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

  rewind(fastq);
  
  string::size_type lpos, rpos;
  while((bytes=getFragment(fastq, &fq))) {
    BOOST_FOREACH(Unmatched& unm, g_unm) {
      lpos = rpos = string::npos;
      lpos = fq.d_nucleotides.find(unm.left);
      if(lpos != string::npos)
	rpos = fq.d_nucleotides.find(unm.right);
      
      if(rpos == string::npos) {
	fq.reverse();
	lpos = fq.d_nucleotides.find(unm.left);
	if(lpos != string::npos)
	  rpos = fq.d_nucleotides.find(unm.right);
      }
      
      if(rpos == string::npos)
	continue;

      // lpos & rpos now correct 
      cout<<"UNM: "<<unm.pos << " - " << unm.pos + unm.unmatched.length() << " (len = "<< unm.unmatched.length()<<", lpos = "<<lpos<<", rpos = "<<rpos<<")"<<endl;
      cout << unm.left << " | "<< unm.unmatched<< " | " << unm.right<<endl;

      Rematched rem;
      if(lpos + 40 > rpos) {
	rem.left = fq.d_nucleotides.substr(lpos, 40);  
	rem.middle="1";
	rem.right=fq.d_nucleotides.substr(rpos);
      }
      else {
	rem.left = fq.d_nucleotides.substr(lpos, 40);
	rem.middle = fq.d_nucleotides.substr(lpos+40, rpos-lpos-40);
	rem.right = fq.d_nucleotides.substr(rpos);
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
}
