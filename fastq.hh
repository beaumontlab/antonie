#pragma once
#include <string>
#include <stdio.h>
#include <stdint.h>
#include <stdexcept>

struct FastQFragment
{
  std::string d_nucleotides;
  std::string d_quality;
  bool exceedsQuality(unsigned int);
  void reverse();
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

  uint64_t getPos() const
  {
    return d_pos;
  }

  void seek(uint64_t pos) 
  {
    fseek(d_fp, pos, SEEK_SET);
    // d_pos gets reset AFTER a read
  }

  unsigned int getFragment(FastQFragment* fq, unsigned int size=0);
private:
  FILE *d_fp;
  uint64_t d_pos;
};

