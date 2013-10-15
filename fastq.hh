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
  bool exceedsQuality(unsigned int);
  void reverse();
  bool reversed;
  uint64_t position;
};

class FASTQReader
{
public:
  FASTQReader(const std::string& str) 
  {
    d_fp = fopen(str.c_str(), "r");
    if(!d_fp) 
      throw std::runtime_error("Unable to open file "+str+" for FASTQ inpuot");
  }

  void seek(uint64_t pos) 
  {
    fseek(d_fp, pos, SEEK_SET);
  }

  unsigned int getRead(FastQRead* fq, unsigned int size=0);
private:
  FILE *d_fp;
};

