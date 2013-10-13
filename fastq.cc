#include "fastq.hh"
#include <stdint.h>
#include <algorithm>
#include "misc.hh"

using namespace std;

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
