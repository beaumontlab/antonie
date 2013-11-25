#pragma once
#include <string>
#include <stdint.h>
#include <vector>
#include <map>

struct Entry16S
{
  uint32_t id;
  int hits;
};

//! A class to match reads to a 16S database
class Search16S
{
public:
  Search16S(const std::string& src, int indexLength);
  void score(const std::string& nucleotides);
  std::vector<Entry16S> topScores();
  void printTop(unsigned int top);
private:
  typedef std::map<uint32_t, std::vector<unsigned int> > hashes_t;
  hashes_t d_hashes;
  std::vector<Entry16S> d_entries;
};
