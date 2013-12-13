#include "dnamisc.hh"
#include "antonie.hh"
#include <vector>
#include <stdexcept>
#include <math.h>
#include <boost/lexical_cast.hpp>

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
    else if(*c=='C') val=2;
    else if(*c=='G') val=2;
    else if(*c=='T') val=3;
    else 
      continue;

    ret |= val;
  }
  return ret;
}
