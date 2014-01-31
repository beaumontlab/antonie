#pragma once
#include <string>
#include <stdint.h>
#include <vector>
#include <map>
#include <stdio.h>
#include "zstuff.hh"

//! A class to match reads to a 16S database
class Search16S
{
public:
  struct Entry
  {
    uint32_t id;
    std::string nucs;
    std::string name;
    bool operator<(const Entry& rhs) const 
    {
      return id < rhs.id;
    }
  };

  Search16S(const std::string& src);
  bool get(Entry* entry);

private:
  std::unique_ptr<LineReader> d_linereader;
};
