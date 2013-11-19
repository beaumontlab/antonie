/* (C) 2013 TU Delft
   (C) 2013 Netherlabs Computer Consulting BV */

#define __STDC_FORMAT_MACROS
#include <tclap/CmdLine.h>
#include <stdio.h>
#include <iostream>
#include <stdint.h>
#include <stdlib.h>
#include <map>
#include <unordered_map>
#include <vector>
#include <string>
#include <string.h>
#include <stdexcept>
#include <forward_list>
#include <inttypes.h>
#include <algorithm>
#include <numeric>
#include <unistd.h>
#include <errno.h>
#include <math.h>
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
#include "geneannotated.hh"
#include "misc.hh"
#include "fastq.hh"

#include <mba/diff.h>
#include <mba/msgno.h>
#include "antonie.hh"
#include "saminfra.hh"

using namespace std;

using namespace boost::algorithm;

extern "C" {
#include "hash.h"
}

namespace io = boost::iostreams;
typedef io::tee_device<std::ostream, std::ostringstream> TeeDevice;
typedef io::stream< TeeDevice > TeeStream;
TeeStream* g_log;

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
  forward_list<FASTQMapping> d_fastqs;
  unsigned int coverage;
};

struct Unmatched
{
  string left, unmatched, right;
  dnapos_t pos;
};


template<typename T>
const typename T::value_type& pickRandom(const T& t)
{
  return t[random() % t.size()];
}

class DuplicateCounter
{
public:
  DuplicateCounter(int estimate=1000000)
  {
    d_hashes.reserve(estimate);
  }
  void feedString(const std::string& str);
  void clear();
  typedef map<uint64_t,uint64_t> counts_t;

  counts_t getCounts();
private:
  vector<uint32_t> d_hashes;
};

void DuplicateCounter::feedString(const std::string& str)
{
  uint32_t hashval = hash(str.c_str(), str.length(), 0);
  d_hashes.push_back(hashval);
}

DuplicateCounter::counts_t DuplicateCounter::getCounts()
{
  counts_t ret;
  sort(d_hashes.begin(), d_hashes.end());
  uint64_t repeatCount=1;
  for(auto iter = next(d_hashes.begin()) ; iter != d_hashes.end(); ++iter) {
    if(*prev(iter) != *iter) {
      ++ret[min(repeatCount, (decltype(repeatCount))20)];
      repeatCount=1;
    }
    else
      repeatCount++;
  }
  ++ret[repeatCount]; 
  return ret;
}

void DuplicateCounter::clear()
{
  d_hashes.clear();
  d_hashes.shrink_to_fit();
}

class ReferenceGenome
{
public:
  ReferenceGenome(const string& fname);
  dnapos_t size() const {
    return d_genome.size() - 1; // we pad at the beginning so we are 1 based..
  }
  vector<uint32_t> getMatchingHashes(const vector<uint32_t>& hashes);
  vector<dnapos_t> getReadPositions(const std::string& nucleotides)
  {
    vector<dnapos_t> ret;
    if(!d_indexes.count(nucleotides.length()))
      throw runtime_error("Attempting to find a read of length we've not indexed for ("+boost::lexical_cast<string>(nucleotides.length())+")");

    auto& index = d_indexes[nucleotides.length()];
  
    uint32_t hashval = hash(nucleotides.c_str(), nucleotides.length(), 0);
    HashPos hp(hashval, 0);
    pair<index_t::const_iterator, index_t::const_iterator> range = equal_range(index.begin(), index.end(), hp);
    if(range.first == range.second)
      return ret;

    for(;range.first != range.second; range.first++) {
      if(!memcmp(d_genome.c_str() + range.first->d_pos, nucleotides.c_str(), nucleotides.length())) {
        ret.push_back(range.first->d_pos);
      }
    }
    return ret;
  }
  
  dnapos_t getReadPosBoth(FastQRead* fq, int qlimit) // tries original & complement
  {
    vector<dnapos_t> positions;

    for(int tries = 0; tries < 2; ++tries) {
      positions = getReadPositions(fq->d_nucleotides);

      if(!positions.empty()) {
	auto pick = pickRandom(positions);

        cover(pick, fq->d_nucleotides.size(), fq->d_quality, qlimit);

        return pick;
      }
      fq->reverse();
    }
    return dnanpos;
  }

  struct MatchDescriptor
  {
    dnapos_t pos;
    bool reverse;
    int score;
  };

  vector<MatchDescriptor> getAllReadPosBoth(FastQRead* fq) // tries original & complement
  {
    vector<MatchDescriptor > ret;
    string nucleotides;
    for(int tries = 0; tries < 2; ++tries) {
      for(auto position : getReadPositions(fq->d_nucleotides)) 
	ret.push_back({position, (bool)tries, 0});
      fq->reverse();
    }
    return ret;
  }

  void cover(dnapos_t pos, unsigned int length, const std::string& quality, int limit) 
  {
    const char* p = quality.c_str();
    for(unsigned int i = 0; i < length; ++i) {
      if(p[i] > limit)
        d_mapping[pos+i].coverage++;
    }
  }

  void cover(dnapos_t pos, char quality, int limit) 
  {
    if(quality > (int) limit)
      d_mapping[pos].coverage++;
  }

  void mapFastQ(dnapos_t pos, const FastQRead& fqfrag, int indel=0)
  {
    FASTQMapping fqm;
    fqm.pos=fqfrag.position;
    fqm.reverse = fqfrag.reversed;
    fqm.indel = indel;
    d_mapping[pos].d_fastqs.push_front(fqm);
    //    cout<<"Adding mapping at pos "<<pos<<", indel = "<<indel<<", reverse= "<<fqm.reverse<<endl;
  }

  vector<dnapos_t> getGCHisto();

  string snippet(dnapos_t start, dnapos_t stop) const { 
    if(stop > d_genome.size()) {
      return d_genome.substr(start);
    }
    return d_genome.substr(start, stop-start);
  }
  void printCoverage(FILE* jsfp, const std::string& fname);
  void index(unsigned int length);

  string getMatchingFastQs(dnapos_t pos, StereoFASTQReader& fastq); 
  string getMatchingFastQs(dnapos_t start, dnapos_t stop,  StereoFASTQReader& fastq); 
  vector<GenomeLocusMapping> d_mapping;
  vector<unsigned int> d_correctMappings, d_wrongMappings, d_gcMappings, d_taMappings;
  vector<vector<uint32_t>> d_kmerMappings;
  vector<Unmatched> d_unmRegions;
  struct LociStats
  {
    vector<std::tuple<char,char,char>> samples; // nucleotide, quality, direction 
  };
  dnapos_t d_aCount, d_cCount, d_gCount, d_tCount;
  typedef unordered_map<dnapos_t, LociStats> locimap_t;
  locimap_t d_locimap;
  unordered_map<dnapos_t, unsigned int> d_insertCounts;
  string d_name;

private:
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
  map<int, index_t> d_indexes;
};

ReferenceGenome::ReferenceGenome(const string& fname)
{
  FILE* fp = fopen(fname.c_str(), "r");
  if(!fp)
    throw runtime_error("Unable to open reference genome file '"+fname+"'");
  d_genome.reserve(filesize(fname.c_str()));  // slight overestimate which is great

  char line[256]="";

  sfgets(line, sizeof(line), fp);
  chomp(line);

  if(line[0] != '>') 
    throw runtime_error("Input not FASTA");
  (*g_log)<<"Reading FASTA reference genome of '"<<line+1<<"'\n";
  char* spacepos=strchr(line+1, ' ');
  if(spacepos)
    *spacepos=0;
  d_name=line+1;

  d_genome="*"; // this gets all our offsets ""right""
  while(fgets(line, sizeof(line), fp)) {
    chomp(line);
    d_genome.append(line);
  }
  
  d_aCount = d_cCount = d_gCount = d_tCount = 0;
  for(auto c : d_genome) {
    acgtDo(c, [&](){ ++d_aCount; }, [&](){ ++d_cCount; }, [&](){ ++d_gCount; }, [&](){ ++d_tCount; });
  }

  d_mapping.resize(d_genome.size());
}

// returns as if we sampled once per index length, an array of index length bins
vector<dnapos_t> ReferenceGenome::getGCHisto()
{
  vector<dnapos_t> ret;
  unsigned int indexlength = d_indexes.rbegin()->first;
  ret.resize(indexlength); // biggest index
  for(dnapos_t pos = 0; pos < d_genome.size() ; pos += indexlength/4) {
    ret[round(indexlength*getGCContent(snippet(pos, pos + indexlength)))]++;
  }
  for(auto& c : ret) {
    c/=4;
  }
  return ret;
}

void ReferenceGenome::index(unsigned int length)
{
  if(length > d_correctMappings.size()) {
    d_correctMappings.resize(length);
    d_wrongMappings.resize(length);
    d_taMappings.resize(length);
    d_gcMappings.resize(length);
    d_kmerMappings.resize(length);
  }

  auto& index = d_indexes[length];
  index.reserve(d_genome.length());
  
  for(string::size_type pos = 0 ; pos < d_genome.length() - length; ++pos) {
    uint32_t hashval = hash(d_genome.c_str() + pos, length, 0);
    index.push_back(HashPos(hashval, pos));
  }

  sort(index.begin(), index.end());
  uint64_t diff = 0;

  for(auto iter = index.begin(); iter!= index.end() ; ++iter) {
    if(iter != index.begin() && iter->d_hash != prev(iter)->d_hash) {
      diff++;
    }
  }
  (*g_log)<<"Average fill in genome hash of length "<<length<<": "<<1.0*d_genome.length()/diff<<endl;
}

string ReferenceGenome::getMatchingFastQs(dnapos_t pos, StereoFASTQReader& fastq)
{
  return getMatchingFastQs(pos > 150 ? pos-150 : 0, pos+150, fastq);
}

string ReferenceGenome::getMatchingFastQs(dnapos_t start, dnapos_t stop, StereoFASTQReader& fastq) 
{
  ostringstream os;

  string reference=snippet(start, stop);
  unsigned int insertPos=0;
  for(unsigned int i = 0 ; i < stop - start; ++i) {
    if(i== (stop-start)/2)
      os << reference << endl;
    string spacer(i, ' ');
    for(auto& fqm : d_mapping[start+i].d_fastqs) {
      FastQRead fqr;
      fastq.getRead(fqm.pos, &fqr);
      if(fqm.reverse)
        fqr.reverse();

      if(fqm.indel > 0 && !insertPos) { // our read has an insert at this position, stretch reference
        if(i+fqm.indel < reference.size())
          reference.insert(i+fqm.indel, 1, '_');
        insertPos=i+fqm.indel;
      } else if(fqm.indel < 0) {      // our read has an erase at this position
        fqr.d_nucleotides.insert(-fqm.indel, 1, 'X');
        fqr.d_quality.insert(-fqm.indel, 1, 42);
      }
      
      if(fqm.indel <= 0 && insertPos && i > insertPos) {
        fqr.d_nucleotides.insert(0, 1, '<');
        fqr.d_quality.insert(0, 1, 40);
      }
      os << spacer;
      int offset=0;
      for(unsigned int j = 0 ; j < fqr.d_nucleotides.size() && i + j + offset < reference.size(); ++j) {
        if(reference[i+j]=='_' && !fqm.indel) {
          os<<'_';
          offset=1;
        }
        if(reference[i+j+offset]==fqr.d_nucleotides[j])
          os<<'.';
        else if(fqr.d_quality[j] > 30) 
          os << fqr.d_nucleotides[j];
        else if(fqr.d_quality[j] < 22) 
          os << ' ';
        else
          os<< (char)tolower(fqr.d_nucleotides[j]);
      }
      os<<"                 "<<(fqm.reverse ? 'R' : ' ');
      os<<endl;
    }
  }
  return os.str();
}

template<typename T>
string jsonVector(const vector<T>& v, const std::string& name, 
		  std::function<double(dnapos_t)> yAdjust = [](dnapos_t d){return 1.0*d;},
		  std::function<double(int)> xAdjust = [](int i){return 1.0*i;})
{
  ostringstream ret;
  ret << "var "<<name<<"=[";
  for(auto iter = v.begin(); iter != v.end(); ++iter) {
    if(iter != v.begin())
      ret<<',';
    ret << '[' << xAdjust(iter - v.begin()) <<','<< yAdjust(*iter)<<']';
  }
  ret <<"];\n";
  return ret.str();
}


template<typename T>
string jsonVectorD(const vector<T>& v, const std::string& name, 
		  std::function<double(double)> yAdjust = [](double d){return d;},
		  std::function<double(int)> xAdjust = [](int i){return 1.0*i;})
{
  ostringstream ret;
  ret << "var "<<name<<"=[";
  for(auto iter = v.begin(); iter != v.end(); ++iter) {
    if(iter != v.begin())
      ret<<',';
    ret << '[' << xAdjust(iter - v.begin()) <<','<< yAdjust(*iter)<<']';
  }
  ret <<"];\n";
  return ret.str();
}

void ReferenceGenome::printCoverage(FILE* jsfp, const std::string& histoName)
{
  uint64_t totCoverage=0, noCoverages=0;
  unsigned int cov;

  bool wasNul=true;
  string::size_type prevNulpos=0;

  vector<unsigned int> covhisto;
  covhisto.resize(1000);
  d_unmRegions.clear();
  for(string::size_type pos = 0; pos < d_mapping.size(); ++pos) {
    cov = d_mapping[pos].coverage;
    bool noCov = cov < 2;
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

  Clusterer<Unmatched> cl(100);
  
  for(auto unm : d_unmRegions) {
    cl.feed(unm);
  }

  (*g_log) << (boost::format("Average depth: %|40t|    %10.2f\n") % (1.0*totCoverage/d_mapping.size())).str();
  (*g_log) << (boost::format("Undercovered nucleotides: %|40t| %10d (%.2f%%), %d ranges\n") % noCoverages % (noCoverages*100.0/d_mapping.size()) % cl.d_clusters.size()).str();

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

uint32_t kmerMapper(const std::string& str, int offset, int unsigned len)
{
  uint32_t ret=0;
  const char *c=str.c_str() + offset;
  string::size_type val;
  for(string::size_type i = 0; i != len; ++i, ++c) {
    ret<<=2;
    if(*c=='A') val=0;
    else if(*c=='C') val=1;
    else if(*c=='G') val=2;
    else if(*c=='T') val=3;
    else 
      continue;

    ret |= val;
  }
  return ret;
}

// 0 if nothing interesting, positive if our read has insert at that position, negative if we have a delete at that position
int MBADiff(dnapos_t pos, const FastQRead& fqr, const string& reference)
{
  string::size_type n, m, d;
  int sn;
  struct varray *ses = varray_new(sizeof(struct diff_edit), NULL);
  
  n = reference.length();
  m = fqr.d_nucleotides.length();
  if ((d = diff(reference.c_str(), 0, n, fqr.d_nucleotides.c_str(), 0, m, NULL, NULL, NULL, 0, ses, &sn, NULL)) == -1) {
    MMNO(errno);
    printf("Error\n");
    return EXIT_FAILURE;
  }
  int ret=0;

#if 0
  if(pos > 3422100 && pos < 3422300) {
    printf("pos %u, d=%lu sn=%d\nUS:  %s\nREF: %s\n", pos, d, 
	   sn, fqr.d_nucleotides.c_str(), reference.c_str());

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
  }
#endif 
  if(sn == 4 && d == 2) {
    struct diff_edit *match1 = (struct diff_edit*)varray_get(ses, 0),
      *change1=(struct diff_edit*)varray_get(ses, 1), 
      *match2=(struct diff_edit*)varray_get(ses, 2), 
      *change2=(struct diff_edit*)varray_get(ses, 3);
    
    if(match1->op == DIFF_MATCH && match2->op==DIFF_MATCH) {
      if(change1->op == DIFF_DELETE && change2->op==DIFF_INSERT) {
	//	cout << "Have delete of "<<change1->len<<" in our read at " << pos+change1->off <<endl;
        ret=-change1->off;
      }
      else if(change1->op == DIFF_INSERT && change2->op==DIFF_DELETE) {
	//        cout<<"Have insert of "<<change1->len<<" in our read at "<<pos+change1->off<<endl;
        ret=change1->off;
      }
    }
  }

  varray_del(ses);

  if(sn > 6)
    return 0;

  return ret;
}

unsigned int diffScore(ReferenceGenome& rg, dnapos_t pos, FastQRead& fqfrag, int qlimit)
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

void printCorrectMappings(FILE* jsfp, const ReferenceGenome& rg, const std::string& name)
{
  fprintf(jsfp, "var %s=[", name.c_str());
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

void printGCMappings(FILE* jsfp, const ReferenceGenome& rg, const std::string& name)
{
  fprintf(jsfp, "var %s=[", name.c_str());
  for(unsigned int i=0; i < rg.d_correctMappings.size() ;++i) {
    double total=rg.d_gcMappings[i] + rg.d_taMappings[i];
    double ratio= rg.d_gcMappings[i]/total;

    fprintf(jsfp, "%s[%d,%.2f]", i ? "," : "", i, ratio);
  }
  fprintf(jsfp,"];\n");
}

struct qtally
{
  qtally() : correct{0}, incorrect{0}{}
  uint64_t correct;
  uint64_t incorrect;
};

int MapToReference(ReferenceGenome& rg, dnapos_t pos, FastQRead fqfrag, int qlimit, SAMWriter* sw, vector<qtally>* qqcounts, int* outIndel=0)
{
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
    if(sw)
      sw->write(pos, fqfrag);
  }
  else {
    int indel=MBADiff(pos, fqfrag, reference);
    if(outIndel)
      *outIndel=indel;
    if(indel) {
      rg.mapFastQ(pos, fqfrag, indel);
      if(sw)
	sw->write(pos, fqfrag, indel);
      didMap=true;
      diffcount=1;
      if(indel > 0) { // our read has an insert at this position
        fqfrag.d_nucleotides.erase(indel, 1); // this makes things align again
        fqfrag.d_quality.erase(indel, 1); 
        rg.d_insertCounts[pos+indel]++;
      } else {      // our read has an erase at this position
        fqfrag.d_nucleotides.insert(-indel, 1, 'X');
        fqfrag.d_quality.insert(-indel, 1, 40);
      }
    }
  }

  //  string diff;
  // diff.reserve(fqfrag.d_nucleotides.length());

  unsigned int readMapPos;
  for(string::size_type i = 0; i < fqfrag.d_nucleotides.size() && i < reference.size();++i) {
    readMapPos = fqfrag.reversed ? ((reference.length()- 1) - i) : i; // d_nucleotides might have an insert
      
    char c =  fqfrag.d_nucleotides[i];

    if(c != reference[i]) {
      //      diff.append(1, fqfrag.d_quality[i] > qlimit ? '!' : '^');
      if(fqfrag.d_quality[i] > qlimit && diffcount < 5) 
        rg.d_locimap[pos+i].samples.push_back(std::make_tuple(fqfrag.d_nucleotides[i], fqfrag.d_quality[i], 
								fqfrag.reversed ^ (i > fqfrag.d_nucleotides.length()/2))); // head or tail
      
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


void emitRegion(FILE*fp, ReferenceGenome& rg, StereoFASTQReader& fastq, GeneAnnotationReader& gar, const string& name, unsigned int index, dnapos_t start, 
        dnapos_t stop, const std::string& report_="")
{
  dnapos_t dnapos = (start+stop)/2;
  fprintf(fp, "region[%d]={name:'%s', pos: %d, depth: [", index, name.c_str(), dnapos);
  for(dnapos_t pos = start; pos < stop; ++pos) {
    if(pos != start) 
      fprintf(fp, ",");
    fprintf(fp, "[%d,%d]", pos, rg.d_mapping[pos].coverage);
  }
  fprintf(fp, "], variability:[");
  for(dnapos_t pos = start; pos < stop; ++pos) {
    if(pos != start) 
      fprintf(fp, ",");
    int count=0;
    if(rg.d_locimap.count(pos)) {
      count=rg.d_locimap[pos].samples.size();
    }
    fprintf(fp, "[%d,%d]", pos, count);
  }

  string picture=rg.getMatchingFastQs(start, stop, fastq);
  replace_all(picture, "\n", "\\n");
  string report = replace_all_copy(report_, "\n", "\\n");

  string annotations;
  auto gas=gar.lookup(dnapos);
  for(auto& ga : gas) {
    annotations += ga.name+" [" + ga.tag  + "], ";
  }
  
  fprintf(fp,"], picture: '%s', annotations: '%s', report: '%s'};\n", picture.c_str(), annotations.c_str(), report.c_str());
  
  fputs("\n", fp);
  fflush(fp);
}

void emitRegion(FILE*fp, ReferenceGenome& rg, StereoFASTQReader& fastq, GeneAnnotationReader& gar, const string& name, unsigned int index, dnapos_t start, const std::string& report="")
{
  emitRegion(fp, rg, fastq, gar, name, index, start-200, start +200, report);
}

unsigned int variabilityCount(const ReferenceGenome& rg, dnapos_t position, const ReferenceGenome::LociStats& lc, double* fraction)
{
  vector<int> counts(256);
  counts[rg.snippet(position, position+1)[0]]+=rg.d_mapping[position].coverage;
  
  int forwardCount=0;

  for(auto& j : lc.samples) {
    counts[get<0>(j)]++;
    if(get<2>(j))
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

vector<ReferenceGenome::MatchDescriptor> fuzzyFind(FastQRead* fqfrag, ReferenceGenome& rg, int keylen, int qlimit)
{
  vector<ReferenceGenome::MatchDescriptor> ret;

  string left, middle, right;
  typedef pair<dnapos_t, char> tpos;
  vector<dnapos_t> lpositions, mpositions, rpositions;

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
			[&match](const ReferenceGenome::MatchDescriptor& md){ return md.pos==match;}) != ret.end())
	  continue;
	
	score = diffScore(rg, match, *fqfrag, qlimit);
        
	ret.push_back({match, fqfrag->reversed, score});
	if(score==0) // won't get any better than this
	  return ret;
      }
    }
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

void printQualities(FILE* jsfp, const qstats_t& qstats)
{
  int i=0;

  fprintf(jsfp, "qualities=[");
  for(const auto& q : qstats) {
    if(i)
      fputs(",", jsfp);
    fprintf(jsfp, "[%d, %f]", i, -10.0*log10(mean(q)));
    ++i;
  }
  fputs("];\n", jsfp);

  vector<double> qlo, qhi;
  for(const auto& q : qstats) {
    qlo.push_back(-10.0*log10(mean(q)) - sqrt(-10.0*log10(variance(q))));
    qhi.push_back(-10.0*log10(mean(q)) +sqrt(-10.0*log10(variance(q))));
  }

  fputs((jsonVectorD(qlo,"qlo")+jsonVectorD(qhi, "qhi")).c_str(), jsfp);

  fflush(jsfp);
}


int main(int argc, char** argv)
{
#ifdef __linux__
  feenableexcept(FE_DIVBYZERO | FE_INVALID); 
#endif 
  
  TCLAP::CmdLine cmd("Command description message", ' ', "g" + string(g_gitHash));

  TCLAP::ValueArg<std::string> annotationsArg("a","annotations","read annotations for reference genome from this file",false, "", "filename", cmd);
  TCLAP::ValueArg<std::string> referenceArg("r","reference","read annotations for reference genome from this file",true,"","string", cmd);
  TCLAP::ValueArg<std::string> fastqArg("f","fastq","read annotations for reference genome from this file",false,"","string", cmd);
  TCLAP::ValueArg<std::string> fastq1Arg("1","fastq1","read annotations for reference genome from this file",true,"","string", cmd);
  TCLAP::ValueArg<std::string> fastq2Arg("2","fastq2","read annotations for reference genome from this file",true,"","string", cmd);
  TCLAP::ValueArg<std::string> excludeArg("x","exclude","read annotations for reference genome from this file",false,"","string", cmd);
  TCLAP::ValueArg<std::string> samFileArg("s","sam-file","Write the assembly to the named SAM file",false,"","filename", cmd);
  TCLAP::ValueArg<int> qualityOffsetArg("q","quality-offset","Quality offset in fastq. 33 for Sanger.",false, 33,"offset", cmd);
  TCLAP::ValueArg<int> beginSnipArg("b","begin-snip","Number of nucleotides to snip from begin of reads",false, 0,"nucleotides", cmd);
  TCLAP::ValueArg<int> endSnipArg("e","end-snip","Number of nucleotides to snip from end of reads",false, 0,"nucleotides", cmd);
  TCLAP::ValueArg<int> qlimitArg("l","qlimit","Disregard nucleotide reads with less quality than this in calls",false, 30,"q", cmd);
  TCLAP::ValueArg<int> duplimitArg("d","duplimit","Ignore reads that occur more than d times. 0 for no filter.",false, 0,"times", cmd);
  TCLAP::SwitchArg unmatchedDumpSwitch("u","unmatched-dump","Create a dump of unmatched reads (unfound.fastq)", cmd, false);

  cmd.parse( argc, argv );

  unsigned int qlimit = qlimitArg.getValue();
  unsigned int duplimit = duplimitArg.getValue();
  srandom(time(0));
  ostringstream jsonlog;  
  TeeDevice td(cerr, jsonlog);
  g_log = new TeeStream(td);

  (*g_log)<<"Antonie was compiled from git hash g" << g_gitHash <<endl;

  GeneAnnotationReader gar(annotationsArg.getValue());
  (*g_log)<<"Done reading "<<gar.size()<<" annotations from '"<<annotationsArg.getValue()<<"'"<<endl;
  
  StereoFASTQReader fastq(fastq1Arg.getValue(), fastq2Arg.getValue(), qualityOffsetArg.getValue(), beginSnipArg.getValue(), endSnipArg.getValue());

  (*g_log)<<"Snipping "<<beginSnipArg.getValue()<<" from beginning of reads, "<<endSnipArg.getValue()<<" from end of reads"<<endl;

  unsigned int bytes=0;
  FastQRead fqfrag1, fqfrag2;
  bytes=fastq.getReadPair(&fqfrag1, &fqfrag2); // get a read to index based on its size
  
  ReferenceGenome rg(referenceArg.getValue());
  double genomeGCRatio = 1.0*(rg.d_cCount + rg.d_gCount)/(rg.d_cCount + rg.d_gCount + rg.d_aCount + rg.d_tCount);
  (*g_log)<<"GC Content of reference genome: "<<100.0*genomeGCRatio<<"%"<<endl;
  rg.index(fqfrag1.d_nucleotides.size());
  
  int keylen=11;
  rg.index(keylen);
  
  unique_ptr<FILE, int(*)(FILE*)> jsfp(fopen("data.js","w"), fclose);

  fprintf(jsfp.get(), "var genomeGCRatio=%f;\n", genomeGCRatio);

  unique_ptr<ReferenceGenome> phix;

  if(!excludeArg.getValue().empty()) {
    phix = unique_ptr<ReferenceGenome>{new ReferenceGenome(excludeArg.getValue())};
    
    (*g_log)<<"Loading positive control filter genome(s)"<<endl;
    
    phix->index(fqfrag1.d_nucleotides.size());
    phix->index(keylen);
  }
  g_log->flush();
  dnapos_t pos;

  uint64_t withAny=0, found=0, total=0, qualityExcluded=0, 
    phixFound=0, differentLength=0, tooFrequent=0, goodPairMatches=0, badPairMatches=0;

  SAMWriter sw(samFileArg.getValue(), rg.d_name, rg.size());

  (*g_log)<<"Performing exact matches of reads to reference genome"<<endl;
  boost::progress_display show_progress(filesize(fastq1Arg.getValue().c_str()), cerr);
 
  for(auto& kmers : rg.d_kmerMappings) 
    kmers.resize(256); // 4^4, corresponds to the '4' below

  qstats_t qstats;
  qstats.resize(fqfrag1.d_nucleotides.size());
  VarMeanEstimator qstat;
  vector<unsigned int> qcounts(256);
  vector<uint64_t> unfoundReads;
  vector<qtally> qqcounts(256);
  vector<dnapos_t> gchisto(fqfrag1.d_nucleotides.size()+1);

  DuplicateCounter dc;
  uint32_t theHash;
  map<uint32_t, uint32_t> seenAlready;
  vector<uint32_t> pairdisthisto(10000);
  do { 
    show_progress += bytes;
    vector<ReferenceGenome::MatchDescriptor > pairpositions[2];
    bool dup1(false), dup2(false);
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
	theHash=hash(fqfrag.d_nucleotides.c_str(), fqfrag.d_nucleotides.size(), 0);
	if(++seenAlready[theHash] > 4) {
	  if(paircount)
	    dup2=true;
	  else
	    dup1=true;
	  tooFrequent++;
	  continue;
	}
      }
      
      gchisto[round(fqfrag.d_nucleotides.size()*getGCContent(fqfrag.d_nucleotides))]++;
      bool hadN=false;
      for(string::size_type i = 0 ; i < fqfrag.d_nucleotides.size(); ++i) {
	char c = fqfrag.d_nucleotides[i];
	if(c=='G' || c=='C')
	  rg.d_gcMappings[i]++;
	else
	  rg.d_taMappings[i]++;
	
	if(fqfrag.d_nucleotides.size() - i > 4)
	  rg.d_kmerMappings[i][kmerMapper(fqfrag.d_nucleotides, i, 4)]++;
	if(c=='N')
	  hadN=true;
      }
      
      if(hadN) {
	unfoundReads.push_back(fqfrag.position);
	withAny++;
	continue;
      }
      /*
      if(fqfrag.d_nucleotides.length() != rg.d_indexlength) {
	differentLength++;
	unfoundReads.push_back(fqfrag.position);
	continue;
      }
      */
      if((pairpositions[paircount]=rg.getAllReadPosBoth(&fqfrag)).empty()) {
	pairpositions[paircount]=fuzzyFind(&fqfrag, rg, keylen, qlimit);
      }
    }
    
    if(pairpositions[0].empty() && pairpositions[1].empty() && phix) {
      auto before=phixFound;
      if(!phix->getAllReadPosBoth(&fqfrag1).empty() || !fuzzyFind(&fqfrag1, *phix, keylen, qlimit).empty())
	phixFound++;
      if(!phix->getAllReadPosBoth(&fqfrag2).empty() || !fuzzyFind(&fqfrag2, *phix, keylen, qlimit).empty())
	phixFound++;
      if(before!=phixFound)
	continue;
    }

    map<int, vector<pair<ReferenceGenome::MatchDescriptor,ReferenceGenome::MatchDescriptor> > > potMatch;
    unsigned int matchCount=0;
    for(auto& match1 : pairpositions[0]) {
      for(auto& match2 : pairpositions[1]) {
	if(match1.reverse != match2.reverse &&  abs((int64_t) match1.pos - (int64_t)match2.pos) < 1400) {
	  potMatch[match1.score + match2.score].push_back({match1, match2});
	  matchCount++;
	}
      }
    }
    if(matchCount > 1) {
      cout<<"Sucks, have "<< matchCount <<" matches for our lovely pair"<<endl;
      for(const auto& scores : potMatch) {
	for(const auto& match : scores.second) {
	  cout<<"\t"<<match.first.pos << " & " << match.second.pos<< " (" << abs((int64_t) match.first.pos - (int64_t)match.second.pos) << ")"<<": "<<match.first.score << " & "<<match.second.score<<endl;
	}
      }
    }
    if(!potMatch.empty()) {
      const auto& chosen = pickRandom(potMatch.begin()->second);
      //	cout<<"Pairwise match!"<<endl;
      goodPairMatches++;
      int distance = chosen.second.reverse ? 
	(fqfrag1.d_nucleotides.length() + (int64_t) chosen.second.pos - (int64_t) chosen.first.pos) :
	(fqfrag1.d_nucleotides.length() + (int64_t) chosen.first.pos - (int64_t) chosen.second.pos);

      if(distance >= 0 && distance < (int)pairdisthisto.size()-1)
	pairdisthisto[distance]++;
      for(int paircount = 0 ; paircount < 2; ++paircount) {
	auto fqfrag = paircount ? &fqfrag2 : &fqfrag1;
	auto dup = paircount ? dup2 : dup1,
	  otherDup = paircount? dup1 : dup2;
	pos = paircount ? chosen.second.pos : chosen.first.pos;

	if((paircount ? chosen.second.reverse : chosen.first.reverse) != fqfrag->reversed)
	  fqfrag->reverse();

	if(otherDup && !dup) {
	  MapToReference(rg, pos, *fqfrag, qlimit, &sw, &qqcounts);
	}
	else if(!otherDup && !dup) {
	  int indel;
	  if(MapToReference(rg, pos, *fqfrag, qlimit, 0, &qqcounts, &indel)) {
	    sw.write(pos, *fqfrag, indel, 3 + (paircount ? 0x80 : 0x40),
		     "=", 
		     paircount ? chosen.first.pos : chosen.second.pos, 
		     (chosen.first.reverse ^ paircount) ? -distance : distance);
	  }
	}
	found++;
      }
    }
    else {
      cout<<"No pair matches, need to map individually: "<<endl;
      badPairMatches++;
      for(unsigned int paircount = 0; paircount < 2; ++paircount) {
	if(paircount ? dup2 : dup1)
	  continue;
	
	map<int, vector<ReferenceGenome::MatchDescriptor>> scores;
	for(auto match: pairpositions[paircount]) {
	  scores[match.score].push_back(match);
	  cout<<"\t"<<paircount<<"\t"<<match.pos<<" "<<match.reverse<<", score: "<<match.score<<endl;
	}
	FastQRead* fqfrag = paircount ? &fqfrag2 : &fqfrag1;
	if(scores.empty()) {
	  unfoundReads.push_back(fqfrag->position);
	  continue;
	}
	auto pick = pickRandom(scores.begin()->second);

	if(fqfrag->reversed != pick.reverse)
	  fqfrag->reverse();

	MapToReference(rg, pick.pos, *fqfrag, qlimit, &sw, &qqcounts);
	found++;
      }
    } 
  } while((bytes=fastq.getReadPair(&fqfrag1, &fqfrag2)));
  
  pairdisthisto.resize(1500);
  fputs(jsonVector(pairdisthisto, "pairdisthisto").c_str(), jsfp.get());

  uint64_t totNucleotides=total*fqfrag1.d_nucleotides.length();
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
  fputs(jsonVector(gchisto, "gcreadhisto",  
		   [totalhisto](dnapos_t c){return 1.0*c/totalhisto;},
		   [&fqfrag1](int i) { return 100.0*i/fqfrag1.d_nucleotides.size();}   ).c_str(), 
	jsfp.get());

  fputs(jsonVector(rg.getGCHisto(), "gcrefhisto",  
		   [&fqfrag1,&rg](dnapos_t c){return 1.0*c/(rg.size()/fqfrag1.d_nucleotides.size());},
		   [&fqfrag1](int i) { return 100.0*i/fqfrag1.d_nucleotides.size();}   ).c_str(), 
	jsfp.get());



  fprintf(jsfp.get(), "var kmerstats=[");
  unsigned int readOffset=0;
  for(const auto& kmer :  rg.d_kmerMappings) {
    if(readOffset >= fqfrag1.d_nucleotides.length() - 4)
      break;

    VarMeanEstimator acc;
    for(auto& count : kmer) {
      acc(count);
    }
    fprintf(jsfp.get(), "%s[%d, %f]", readOffset ? "," : "", readOffset, sqrt(variance(acc)) / mean(acc) );
    readOffset++;
  }
  fprintf(jsfp.get(), "];\n");

  printGCMappings(jsfp.get(), rg, "gcRatios");

  (*g_log) << (boost::format("Total reads: %|40t| %10d (%.2f gigabps)") % total % (totNucleotides/1000000000.0)).str() <<endl;
  (*g_log) << (boost::format("Excluded control reads: %|40t|-%10d") % phixFound).str() <<endl;
  (*g_log) << (boost::format("Quality excluded: %|40t|-%10d") % qualityExcluded).str() <<endl;
  (*g_log) << (boost::format("Ignored reads with N: %|40t|-%10d") % withAny).str()<<endl;
  if(duplimit)
    (*g_log) << (boost::format("Too frequent reads: %|40t| %10d (%.02f%%)") % tooFrequent % (100.0*tooFrequent/total)).str() <<endl;
  (*g_log) << (boost::format("Different length reads: %|40t|-%10d") % differentLength).str() <<endl;
  (*g_log) << (boost::format("Full matches: %|40t|-%10d (%.02f%%)\n") % found % (100.0*found/total)).str();
  (*g_log) << (boost::format(" Reads matched in a good pair: %|40t| %10d\n") % (goodPairMatches*2)).str();
  (*g_log) << (boost::format(" Reads not matched, bad pair: %|40t| %10d\n") % (badPairMatches*2)).str();

  (*g_log) << (boost::format("Not fully matched: %|40t|=%10d (%.02f%%)\n") % unfoundReads.size() % (unfoundReads.size()*100.0/total)).str();
  (*g_log) << (boost::format("Mean Q: %|40t|    %10.2f +- %.2f\n") % (-10.0*log10(mean(qstat))) 
	       % sqrt(-10.0*log10(variance(qstat)) )).str();

  seenAlready.clear();

  for(auto& i : rg.d_correctMappings) {
    i=found;
  }

  if(phix) {
    for(auto& i : phix->d_correctMappings) {
      i=phixFound;
    }
  }

  rg.printCoverage(jsfp.get(), "fullHisto");
  printQualities(jsfp.get(), qstats);

  if(unmatchedDumpSwitch.getValue())
    writeUnmatchedReads(unfoundReads, fastq);

  int index=0;

  Clusterer<Unmatched> cl(100);
  for(auto unm : rg.d_unmRegions) {
    cl.feed(unm);
  }

  for(auto unmCl : cl.d_clusters) {
    emitRegion(jsfp.get(), rg, fastq, gar, "Undermatched", index++, unmCl.getBegin()-100, unmCl.getEnd()+100);
  }
  
  printCorrectMappings(jsfp.get(), rg, "referenceQ");
  if(phix)
    printCorrectMappings(jsfp.get(), *phix, "controlQ");
  fprintf(jsfp.get(), "var qqdata=[");
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

      fprintf(jsfp.get(), "%s[%ld, %f, %" PRIu64 "]", 
	      printedYet ?  "," : "", 
	      coinco - qqcounts.begin(), qscore,
	      coinco->incorrect + coinco->correct);
      printedYet=true;
    }
  }
  fprintf(jsfp.get(), "];\n");


  uint64_t significantlyVariable=0;
  boost::format fmt1("%-10d: %3d*%c ");
  string fmt2("                  ");
  int aCount, cCount, tCount, gCount;
  double fraction;
  
  struct ClusterLocus
  {
    unsigned int pos;
    ReferenceGenome::LociStats locistat;
    bool operator<(const ClusterLocus& rhs) const
    {
      return pos < rhs.pos;
    }
  };
  vector<ClusterLocus> clusters;

  for(auto& locus : rg.d_locimap) {
    clusters.push_back(ClusterLocus{locus.first, locus.second});
  }

  sort(clusters.begin(), clusters.end());

  Clusterer<ClusterLocus> cll(100);
  for(auto item : clusters) 
    cll.feed(item);

  (*g_log)<<"Found "<<rg.d_locimap.size()<<" varying loci ("<<cll.d_clusters.size()<<" ranges)"<<endl;  

  for(auto& cluster : cll.d_clusters) {
    vector<pair<string, dnapos_t>> reports;
    for(auto& locus : cluster.d_members) {
      unsigned int varcount=variabilityCount(rg, locus.pos, locus.locistat, &fraction);
      if(varcount < 10) 
	continue;
      ostringstream report;

      char c=rg.snippet(locus.pos, locus.pos+1)[0];
      aCount = cCount = tCount = gCount = 0;

      acgtDo(c, 
	     [&](){aCount += rg.d_mapping[locus.pos].coverage;},
	     [&](){cCount += rg.d_mapping[locus.pos].coverage;},
	     [&](){gCount += rg.d_mapping[locus.pos].coverage;},
	     [&](){tCount += rg.d_mapping[locus.pos].coverage;}
	     );
      report << (fmt1 % locus.pos % rg.d_mapping[locus.pos].coverage % rg.snippet(locus.pos, locus.pos+1) ).str();
      sort(locus.locistat.samples.begin(), locus.locistat.samples.end());

      significantlyVariable++;
      for(auto j = locus.locistat.samples.begin(); 
	  j != locus.locistat.samples.end(); ++j) {
	c=get<0>(*j);
	acgtDo(c, [&](){ aCount++; }, [&](){ cCount++; }, [&](){ gCount++; }, [&](){ tCount++; } );
      
	report<<c;
      }
      report<<endl<<fmt2;
      for(auto j = locus.locistat.samples.begin(); 
	  j != locus.locistat.samples.end(); ++j) {
	report<<((char)(get<1>(*j)+33));
      }
      report << endl << fmt2;
      for(auto j = locus.locistat.samples.begin(); 
	  j != locus.locistat.samples.end(); ++j) {
	report << (get<2>(*j) ? 'R' : '.');
      }

      int tot=locus.locistat.samples.size() + rg.d_mapping[locus.pos].coverage;
      report<<endl;
      vector<GeneAnnotation> gas=gar.lookup(locus.pos);
      if(!gas.empty()) {
	report << fmt2 << "Annotation: ";
	for(auto& ga : gas) {
	  report << ga.name<<" ["<<ga.tag<<"], ";
	}
	report << endl;
      }
      report << fmt2<< "Fraction tail: "<<fraction<<", "<< locus.locistat.samples.size()<<endl;
      report << fmt2<< "A: " << aCount*100/tot <<"%, C: "<<cCount*100/tot<<"%, G: "<<gCount*100/tot<<"%, T: "<<tCount*100/tot<<"%"<<endl;

      // cout<<rg.getMatchingFastQs(locus.pos, fastq);
      reports.push_back({report.str(), locus.pos});  
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

      emitRegion(jsfp.get(), rg, fastq, gar, "Variable", index++, minPos-100, maxPos+100, theReport);
    }
  }
  (*g_log)<<"Found "<<significantlyVariable<<" significantly variable loci"<<endl;
  (*g_log)<<"Found "<<rg.d_insertCounts.size()<<" loci with at least one insert in a read"<<endl;
  struct revsort
  {
    bool operator()(const unsigned int&a, const unsigned int&b) const
    { return a > b;} 
  };
  map<unsigned int, vector<dnapos_t>, revsort> topInserts;
  
  unsigned int significantInserts=0;
  for(const auto& insloc : rg.d_insertCounts) {
    topInserts[insloc.second].push_back(insloc.first);
    if(insloc.second > 10)
      significantInserts++;
  }
  (*g_log)<<"Found "<<significantInserts<<" significant inserts"<<endl;
  
  for(const auto& insert : topInserts) {
    if(insert.first < 10)
      break;
    for(const auto& position : insert.second) {
      cout<<position<<"\t"<<insert.first<<" inserts"<<endl;
      vector<GeneAnnotation> gas=gar.lookup(position);
      if(!gas.empty()) {
	cout<<fmt2<<"Annotation: ";
	for(auto& ga : gas) {
	  cout<<ga.name<<" ["<<ga.tag<<"], ";
	}
	cout<<endl;
      }
      emitRegion(jsfp.get(), rg, fastq, gar, "Insert", index++, position);
      cout<<rg.getMatchingFastQs(position, fastq);
    }
  }
  
  g_log->flush();
  string log = jsonlog.str();
  replace_all(log, "\n", "\\n");
  fprintf(jsfp.get(), "var antonieLog=\"%s\";\n", log.c_str());

  exit(EXIT_SUCCESS);
}

#if 0
  /*
  rg.index(75);
  auto chimericFound=halfFind(unfoundReads, fastq, rg, 75);
  (*g_log)<<(boost::format("Probable chimeric reads:%|40t|-%10d (%.2f%%)\n") 
         % chimericFound % (100.0*chimericFound/total)).str();
  */

  /*
  rg.printFastQs(5718000, fastq);  
  */
  /*
  */
#endif
