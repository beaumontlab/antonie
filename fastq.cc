#include "fastq.hh"
#include <stdint.h>
#include <algorithm>
#include <string.h>
#include "misc.hh"

using namespace std;

FASTQReader::FASTQReader(const std::string& str) 
{
  d_fp = fopen(str.c_str(), "r");
  if(!d_fp) 
    throw std::runtime_error("Unable to open file '"+str+"' for FASTQ input");
}

bool FastQRead::exceedsQuality(unsigned int limit)
{
  uint8_t q;
  for(string::size_type pos = 0 ; pos < d_quality.size(); ++pos) {
    q = d_quality[pos] - 33;
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

unsigned int FASTQReader::getRead(FastQRead* fq, unsigned int size)
{
  uint64_t pos = ftell(d_fp);
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
  fq->d_quality=line;
  fq->reversed=0;
  fq->position=pos;
  return ftell(d_fp) - pos;
}
