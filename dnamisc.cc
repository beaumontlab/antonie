#include "dnamisc.hh"
#include "antonie.hh"
#include <vector>
#include <stdexcept>
#include <math.h>
extern "C" {
#include "hash.h"
}
#include <algorithm>
#include <boost/lexical_cast.hpp>
dnapos_t dnanpos = (dnapos_t) -1;
using std::vector;
using std::runtime_error;

double getGCContent(const std::string& str)
{
  dnapos_t aCount{0}, cCount{0}, gCount{0}, tCount{0}, nCount{0};
  for(auto c : str) {
    if(c=='A') ++aCount;
    else if(c=='C') ++cCount;
    else if(c=='G') ++gCount;
    else if(c=='T') ++tCount;
    else if(c=='N') ++nCount;
  }
  dnapos_t total = cCount + gCount + aCount + tCount + nCount;
  if(!total)
    return 0.0;
  return 1.0*(cCount + gCount)/(1.0*total);
}

double qToErr(unsigned int i) 
{
  static vector<double> answers;
  
  if(answers.empty()) {
    for(int n = 0; n < 60 ; ++n) {
      answers.push_back(pow(10.0, -n/10.0));
    }
  }
  if(i > answers.size()) {
    throw runtime_error("Can't calculate error rate for Q "+boost::lexical_cast<std::string>(i));
  }

  return answers[i];
}

uint32_t kmerMapper(const std::string& str, int offset, int unsigned len)
{
  uint32_t ret=0;
  const char *c=str.c_str() + offset;
  std::string::size_type val;
  for(std::string::size_type i = 0; i != len; ++i, ++c) {
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

const char* AminoAcidName(char c)
{
  switch(c) {
    case 'T':
      return "Threonine";
    case 'F':
      return "Phenylanaline";
    case 'L':
      return "Leucine";
    case 'I':
      return "Isoleucine";
    case 'M':
      return "Methionine";
    case 'V':
      return "Valine";
    case 'S':
      return "Serine";
    case 'P':
      return "Proline";
    case 'A':
      return "Alanine";
    case 'Y':
      return "Tyrosine";
    case 's':
      return "Stop";
    case 'H':
      return "Histidine";
    case 'Q':
      return "Glutamine";
    case 'N':
      return "Asparagine";
    case 'K':
      return "Lysine";
    case 'D':
      return "Aspartic Acid";
    case 'E':
      return "Glutamic Acid";
    case 'C':
      return "Cysteine";
    case 'W':
      return "Tryptophan";
    case 'R':
      return "Arganine";
    case 'G':
      return "Glycine";
  }
  return "?";
  
}

char DNAToAminoAcid(const char* s)
{
  char a=*s++;
  char b=*s++;
  char c=*s;
  if(a=='T') {
    if(b=='T')  {
      if(c=='T' || c=='C')
        return 'F';
      else 
        return 'L';
    }
    if(b=='C')
      return 'S';
    if(b=='A') {
      if(c=='T' || c=='C')
        return 'Y';
      else
        return 's';
    }
    if(b=='G') {
      if(c=='T' || c=='C')
        return 'C';
      else if(c=='A')
        return 's';
      else if(c=='G')
        return 'W';
    }
  }
  else if(a=='C') {
    if(b=='T')
      return 'L';
    if(b=='C')
      return 'P';
    if(b=='A') {
      if (c=='T' || c=='C')
        return 'H';
      else  
        return 'Q';
    }
    if(b=='G')
      return 'R';
  }
  else if(a=='A') {
    if(b=='T') {
      if(c=='G')
        return 'M';
      else
        return 'I';
    }
    if(b=='C')
      return 'T';
    else if(b=='A') {
      if(c=='T' || c=='C')
        return 'N';
      else
        return 'K';
    }
    else if(b=='G') {
      if(c=='T' || c=='C')
        return 'S';
      else
        return 'R';
    }
  }
  else if(a=='G') {
    if(b=='T')
      return 'V';
    else if(b=='C')
      return 'A';
    else if(b=='A') {
      if(c=='T' || c=='C')
        return 'D';
      else
        return 'E';
    }
    else if(b=='G')
      return 'G';
  }
  return '?';
}

void DuplicateCounter::feedString(const std::string& str)
{
  uint32_t hashval = qhash(str.c_str(), str.length(), 0);
  d_hashes.push_back(hashval);
}

DuplicateCounter::counts_t DuplicateCounter::getCounts()
{
  counts_t ret;
  sort(d_hashes.begin(), d_hashes.end());
  uint64_t repeatCount=1;
  for(auto iter = next(d_hashes.begin()) ; iter != d_hashes.end(); ++iter) {
    if(*prev(iter) != *iter) {
      ret[std::min(repeatCount, (decltype(repeatCount))20)]+=repeatCount; 
      repeatCount=1;
    }
    else
      repeatCount++;
  }
  ret[repeatCount]+=repeatCount;
  return ret;
}

void DuplicateCounter::clear()
{
  d_hashes.clear();
  d_hashes.shrink_to_fit();
}
