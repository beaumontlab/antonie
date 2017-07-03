#include "refgenome.hh"
#include <boost/lexical_cast.hpp>
#include <stdexcept>
#include <string.h>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include "misc.hh"
#include "dnamisc.hh"

extern "C" {
#include "hash.h"
}

using boost::lexical_cast;
using namespace std;

vector<dnapos_t> ReferenceChromosome::getReadPositions(const std::string& nucleotides)
{
  vector<dnapos_t> ret;
  if(!d_indexes.count(nucleotides.length()))
    throw runtime_error("Attempting to find a read of length we've not indexed for ("+boost::lexical_cast<string>(nucleotides.length())+")");
  
  auto& index = d_indexes[nucleotides.length()];
  
  uint32_t hashval = qhash(nucleotides.c_str(), nucleotides.length(), 0);
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
  
dnapos_t ReferenceChromosome::getReadPosBoth(FastQRead* fq, int qlimit) // tries original & complement
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

vector<ReferenceChromosome::MatchDescriptor> ReferenceChromosome::getAllReadPosBoth(FastQRead* fq) // tries original & complement
{
  vector<MatchDescriptor > ret;
  string nucleotides;
  for(int tries = 0; tries < 2; ++tries) {
    for(auto position : getReadPositions(fq->d_nucleotides)) 
      ret.push_back({this, position, (bool)tries, 0});
    fq->reverse();
  }
  return ret;
}

void ReferenceChromosome::cover(dnapos_t pos, unsigned int length, const std::string& quality, int limit) 
{
  const char* p = quality.c_str();
  for(unsigned int i = 0; i < length; ++i) {
    if(p[i] > limit)
      d_mapping[pos+i].coverage++;
  }
}

void ReferenceChromosome::cover(dnapos_t pos, char quality, int limit) 
{
  if(quality > (int) limit)
    d_mapping[pos].coverage++;
}

void ReferenceChromosome::mapFastQ(dnapos_t pos, const FastQRead& fqfrag, int indel)
{
  FASTQMapping fqm;
  fqm.pos=fqfrag.position;
  fqm.reverse = fqfrag.reversed;
  fqm.indel = indel;
  d_mapping[pos].d_fastqs.push_front(fqm);
  //    cout<<"Adding mapping at pos "<<pos<<", indel = "<<indel<<", reverse= "<<fqm.reverse<<endl;
}



string ReferenceChromosome::snippet(dnapos_t start, dnapos_t stop) const 
{ 
  if(stop > d_genome.size()) {
      return d_genome.substr(start);
  }
  return d_genome.substr(start, stop-start);
}

ReferenceChromosome::ReferenceChromosome(const string& fname)
{
  FILE* fp = fopen(fname.c_str(), "rb");
  if(!fp)
    throw runtime_error("Unable to open reference genome file '"+fname+"'");
  d_genome.reserve(filesize(fname.c_str()));  // slight overestimate which is great

  char line[256]="";

  sfgets(line, sizeof(line), fp);
  chomp(line);

  if(line[0] != '>') 
    throw runtime_error("Input not FASTA");
  d_fullname=line+1;

  char* spacepos=strchr(line+1, ' ');

  if(spacepos)
    *spacepos=0;
  d_name=line+1;
  
  d_genome="*"; // this gets all our offsets ""right""
  while(fgets(line, sizeof(line), fp)) {
    chomp(line);
    d_genome.append(line);
  }
  
  initGenome();
}

unique_ptr<ReferenceChromosome> ReferenceChromosome::makeFromString(const std::string& genome)
{
  istringstream istr(genome);
  unique_ptr<ReferenceChromosome> ret(new ReferenceChromosome);
  getline(istr, ret->d_name);
  if(ret->d_name.empty() || ret->d_name[0]!='>') 
    throw runtime_error("Input not FASTA");
  ret->d_fullname=ret->d_name.substr(1); // skip >

  ret->d_name = ret->d_fullname; // should stop after ' '
  auto spacepos = ret->d_name.find(' ');
  if(spacepos != string::npos)
    ret->d_name=ret->d_name.substr(0, spacepos);

  string line;
  ret->d_genome="*"; // this gets all our offsets ""right""
  while(getline(istr, line)) {
    boost::trim_right(line);
    ret->d_genome.append(line);
  }
  
  ret->initGenome();
  return ret;
}

void ReferenceChromosome::initGenome()
{
  d_aCount = d_cCount = d_gCount = d_tCount = 0;
  for(const auto& c : d_genome) {
    switch(c) {
    case 'A':
      ++d_aCount;
      break;
    case 'C':
      ++d_cCount;
      break;
    case 'G':
      ++d_gCount;
      break;
    case 'T':
      ++d_tCount;
      break;
    }
  }

  d_mapping.resize(d_genome.size());
}

// returns as if we sampled once per index length, an array of index length bins
vector<dnapos_t> ReferenceChromosome::getGCHisto()
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

void ReferenceChromosome::index(unsigned int length)
{
  if(length > d_correctMappings.size()) {
    d_correctMappings.resize(length);
    d_wrongMappings.resize(length);
    d_taMappings.resize(length);
    d_gcMappings.resize(length);
  }

  auto& index = d_indexes[length];
  index.reserve(d_genome.length());
  
  for(string::size_type pos = 0 ; pos < d_genome.length() - length; ++pos) {
    uint32_t hashval = qhash(d_genome.c_str() + pos, length, 0);
    index.push_back(HashPos(hashval, pos));
  }

  sort(index.begin(), index.end());
  uint64_t diff = 0;

  for(auto iter = index.begin(); iter!= index.end() ; ++iter) {
    if(iter != index.begin() && iter->d_hash != prev(iter)->d_hash) {
      diff++;
    }
  }
  // (*g_log)<<"Average fill in genome hash of length "<<length<<": "<<1.0*d_genome.length()/diff<<endl;
}

string ReferenceChromosome::getMatchingFastQs(dnapos_t pos, StereoFASTQReader& fastq)
{
  return getMatchingFastQs(pos > 150 ? pos-150 : 1, pos+150, fastq);
}

string ReferenceChromosome::getMatchingFastQs(dnapos_t start, dnapos_t stop, StereoFASTQReader& fastq) 
{
  ostringstream os;
  if(stop > size())
    stop = size();
  if(start > size())
    start = 1;
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

