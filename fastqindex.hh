#pragma once
#include <vector>
#include <memory>
#include "fastq.hh"


struct HashedPos
{
  uint32_t hash;
  uint64_t position;
  bool operator<(const HashedPos& b) const {
    return hash < b.hash;
  } 
  bool operator<(uint32_t h) const {
    return hash < h;
  }
}__attribute__((packed));

std::unique_ptr<std::vector<HashedPos> > indexFASTQ(FASTQReader* fqreader, const std::string& fname, int chunklen);

std::vector<FastQRead> getConsensusMatches(const std::string& consensus, const std::map<FASTQReader*, std::unique_ptr<std::vector<HashedPos> > >& fhpos, int chunklen);
