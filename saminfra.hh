#pragma once
#include <string>
#include <vector>
#include "antonie.hh"
#include <stdio.h>
#include "fastq.hh"
#include "zstuff.hh"

//! Write SAM files, with support for paired-end read mappings
class SAMWriter
{
public:
  SAMWriter(const std::string& fname, const std::string& genome, dnapos_t len);
  ~SAMWriter();
  void write(dnapos_t pos, const FastQRead& fqfrag, int indel=0, int flags=0, const std::string& rnext="*", dnapos_t pnext=0, int32_t tlen=0 );
private:
  FILE* d_fp;
  std::string d_fname;
  std::string d_genomeName;
};


//! Write BAM files, with support for paired-end read mappings
class BAMWriter
{
public:
  BAMWriter(const std::string& fname, const std::string& genome, dnapos_t len);
  ~BAMWriter();
  uint64_t write(dnapos_t pos, const FastQRead& fqfrag, int indel=0, int flags=0, const std::string& rnext="*", dnapos_t pnext=0, int32_t tlen=0 );
  void qwrite(dnapos_t pos, const FastQRead& fqfrag, int indel=0, int flags=0, const std::string& rnext="*", dnapos_t pnext=0, int32_t tlen=0 );
  void runQueue(StereoFASTQReader& sfq);
private:

  std::string d_fname;
  std::string d_genomeName;
  BGZFWriter d_zw;
  FILE* d_baifp;
  struct Write
  {
    bool operator<(const Write& rhs) const
    {
      return pos < rhs.pos;
    }
    dnapos_t pos;
    uint64_t fpos;
    bool reversed;
    int indel;
    int flags;
    std::string rnext;
    dnapos_t pnext;
    int tlen;
    uint64_t voffset;
    unsigned int bin;
  };
  std::vector<Write> d_queue;
};

std::string bamCompress(const std::string& dna);
