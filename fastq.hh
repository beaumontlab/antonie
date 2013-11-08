#pragma once
#include <string>
#include <stdio.h>
#include <stdint.h>
#include <stdexcept>

struct FastQRead
{
  FastQRead() : reversed(false), position(0) {}
  std::string d_nucleotides;
  std::string d_quality;
  std::string d_header;
  bool exceedsQuality(unsigned int);
  std::string getSangerQualityString() const;
  void reverse();
  bool reversed;
  uint64_t position;
};

class FASTQReader
{
public:
  FASTQReader(const std::string& str, unsigned int qoffset, unsigned int snipLeft=0, unsigned int snipRight=0);

  void seek(uint64_t pos) 
  {
    fseek(d_fp, pos, SEEK_SET);
  }

  unsigned int getRead(FastQRead* fq);
private:
  FILE *d_fp;
  unsigned int d_qoffset;
  unsigned int d_snipLeft, d_snipRight;
};

