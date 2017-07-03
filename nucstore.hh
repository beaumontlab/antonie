#pragma once
#include <string>
#include <boost/utility/string_ref.hpp>
extern "C" {
#include "hash.h"
}

class NucleotideStore
{
public:
  void append(char c);
  void append(const boost::string_ref& line);
  char get(size_t pos) const;
  NucleotideStore getRange(size_t pos, size_t len) const;
  size_t size() const
  {
    return 4*d_storage.size() + bitpos/2;
  }

  size_t hash() const
  {
    return qhash(d_storage.c_str(), d_storage.size(), bitpos ? d_curval : 0);
  }

  bool operator==(const NucleotideStore& rhs) const
  {
    return d_storage == rhs.d_storage && bitpos == rhs.bitpos && d_curval == rhs.d_curval;
  }
private:
  int bitpos{0};
  std::string d_storage;
  uint8_t d_curval{0};
};

std::ostream& operator<<(std::ostream& os, const NucleotideStore& ns);
