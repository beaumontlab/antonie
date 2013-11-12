#include "dnamisc.hh"
#include "antonie.hh"

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
