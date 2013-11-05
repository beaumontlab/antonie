#pragma once
#include <string>
#include "antonie.hh"
#include <stdio.h>
#include "fastq.hh"

class SAMWriter
{
public:
  SAMWriter(const std::string& fname, const std::string& genome, dnapos_t len);
  ~SAMWriter();
  void write(dnapos_t dnapos, const FastQRead& fqfrag, int indel=0);

private:
  FILE* d_fp;
  std::string d_fname;
  std::string d_genomeName;
};
