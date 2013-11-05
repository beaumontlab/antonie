#include "fastq.hh"
#include <stdint.h>
#include <algorithm>
#include <string.h>
#include "misc.hh"
#include <boost/lexical_cast.hpp>
using namespace std;

FASTQReader::FASTQReader(const std::string& str, unsigned int qoffset, unsigned int snipLeft, unsigned int snipRight) 
  :  d_snipLeft{snipLeft}, d_snipRight{snipRight} 
{
  d_fp = fopen(str.c_str(), "r");
  if(!d_fp) 
    throw std::runtime_error("Unable to open file '"+str+"' for FASTQ input");
  d_qoffset=qoffset;
}

bool FastQRead::exceedsQuality(unsigned int limit)
{
  uint8_t q;
  for(string::size_type pos = 0 ; pos < d_quality.size(); ++pos) {
    q = d_quality[pos];
    if(q < limit)
      return false;
  }
  return true;
}

void FastQRead::reverse()
{
  reverseNucleotides(&d_nucleotides);
  std::reverse(d_quality.begin(), d_quality.end());
  reversed = !reversed;
}

unsigned int FASTQReader::getRead(FastQRead* fq)
{
  uint64_t pos = ftell(d_fp);
  char line[256]="";
  if(!fgets(line, sizeof(line), d_fp)) 
    return 0;
  if(line[0] != '@')
    throw runtime_error("Input not FASTQ");

  chomp(line);
  fq->d_header.assign(line+1);

  sfgets(line, sizeof(line), d_fp);
  chomp(line);
  
  if(d_snipLeft || d_snipRight)
    fq->d_nucleotides.assign(line + d_snipLeft, strlen(line) -d_snipLeft-d_snipRight);
  else
    fq->d_nucleotides.assign(line);
  sfgets(line, sizeof(line), d_fp);
  sfgets(line, sizeof(line), d_fp);
  chomp(line);

  if(d_snipLeft || d_snipRight)
    fq->d_quality.assign(line + d_snipLeft, strlen(line)-d_snipLeft-d_snipRight);
  else
    fq->d_quality.assign(line);

  for(auto& c : fq->d_quality) {
    if((unsigned int)c < d_qoffset)
      throw runtime_error("Attempting to parse a quality code of val "+boost::lexical_cast<string>((int)c)+" which is < our quality offset");
    c -= d_qoffset;
  }

  fq->reversed=0;
  fq->position=pos;
  return ftell(d_fp) - pos;
}
