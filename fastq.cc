#include "fastq.hh"
#include <stdint.h>
#include <algorithm>
#include <string.h>
#include "misc.hh"
#include <boost/lexical_cast.hpp>
using namespace std;

constexpr uint64_t StereoFASTQReader::s_mask;

FASTQReader::FASTQReader(const std::string& str, unsigned int qoffset, unsigned int snipLeft, unsigned int snipRight) 
  :  d_snipLeft{snipLeft}, d_snipRight{snipRight}, d_reader(LineReader::make(str)) 
{
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

string FastQRead::getSangerQualityString() const
{
  string quality{d_quality};
  for(auto& c : quality) 
    c+=33;
  return quality;
}

void FastQRead::reverse()
{
  reverseNucleotides(&d_nucleotides);
  std::reverse(d_quality.begin(), d_quality.end());
  reversed = !reversed;
}

unsigned int FASTQReader::getRead(FastQRead* fq)
{
  uint64_t pos = d_reader->getUncPos();
  char line[1024]="";
  if(!d_reader->fgets(line, sizeof(line)))
    return 0;
  if(line[0] != '@')
    throw runtime_error("Input not FASTQ");

  chomp(line);
  fq->d_header.assign(line+1);

  d_reader->fgets(line, sizeof(line));
  chomp(line);
  
  if(d_snipLeft || d_snipRight)
    fq->d_nucleotides.assign(line + d_snipLeft, strlen(line) -d_snipLeft-d_snipRight);
  else
    fq->d_nucleotides.assign(line);
  d_reader->fgets(line, sizeof(line));
  d_reader->fgets(line, sizeof(line));

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
  return d_reader->getUncPos() - pos;
}

unsigned int StereoFASTQReader::getRead(uint64_t pos, FastQRead* fq)
{
  unsigned int ret;
  if(pos & (1ULL<<63)) {
    d_fq2.seek(pos & s_mask);
    ret=d_fq2.getRead(fq);
  }
  else {
    d_fq1.seek(pos);
    ret=d_fq1.getRead(fq);
  }

  fq->position = pos;
  return ret;
}

unsigned int StereoFASTQReader::getReadPair(FastQRead* fq1, FastQRead* fq2)
{
  unsigned int ret1, ret2;
  ret1=d_fq1.getRead(fq1);
  ret2=d_fq2.getRead(fq2);

  if(ret1 != ret2) {
    throw runtime_error("Difference between paired files in read: " + boost::lexical_cast<string>(ret1) +" != "+ boost::lexical_cast<string>(ret2));
  }
  fq2->position |= (1ULL<<63);
  return ret1;
}



