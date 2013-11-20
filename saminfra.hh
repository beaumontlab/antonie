#pragma once
#include <string>
#include "antonie.hh"
#include <stdio.h>
#include <vector>
#include <algorithm>
#include "fastq.hh"

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

class SAMSorter
{
public:
  void put(uint64_t foffset, dnapos_t pos)
  {
    d_entries.push_back({foffset,pos});
  }
  std::vector<uint64_t> getSortedFOffsets()
  {
    sort(d_entries.begin(), d_entries.end());
    std::vector<uint64_t> ret;
    ret.reserve(d_entries.size());
    for(auto e: d_entries) {
      ret.push_back(e.foffset);
    }
    return ret;
  }
private:
  struct entry {
    uint64_t foffset;
    dnapos_t pos;
    bool operator<(const entry& b) const
    {
      return (pos & (~(1ULL<<63))) < (b.pos & (~(1ULL<<63)));
    }
  };
  std::vector<entry> d_entries;
};
