#pragma once
#include <string>
#include <vector>
#include <boost/utility/string_ref.hpp>
#include <string.h>
extern "C" {
#include "hash.h"
}

class NucleotideStore
{
public:
  explicit NucleotideStore(const boost::string_ref& in)
  {
    append(in);
  }
  NucleotideStore() {}
  void append(char c);
  void append(const boost::string_ref& line);
  char get(size_t pos) const;
  char operator[](size_t pos) const
  {
    return get(pos);
  }
  void set(size_t pos, char c);
  NucleotideStore getRange(size_t pos, size_t len) const;
  NucleotideStore getRC() const;
  size_t size() const
  {
    return 4*d_storage.size() + bitpos/2;
  }

  struct Delta
  {
    uint32_t pos;
    char o;
    enum class Action {Replace, Delete, Insert} a;
    bool operator==(const Delta& rhs) const
    {
      return pos==rhs.pos && o==rhs.o && a==rhs.a;
    }
  };

  std::vector<Delta> getDelta(const NucleotideStore& b, double mispen=1, double gappen=2, double skwpen=0) const;
  void applyDelta(std::vector<Delta>& delta);
  size_t hash() const
  {
    /*
    if(d_storage.size()==4) {
      uint32_t ret;
      memcpy((char*)&ret, d_storage.c_str(), 4);
      return ret;
    }
    */
    return qhash(d_storage.c_str(), d_storage.size(), bitpos ? d_curval : 0);
  }

  size_t overlap(const NucleotideStore& rhs) const;
  size_t fuzOverlap(const NucleotideStore& rhs, int ratio) const;
  
  bool operator==(const NucleotideStore& rhs) const
  {
    return d_storage == rhs.d_storage && bitpos == rhs.bitpos && d_curval == rhs.d_curval;
  }

  bool operator<(const NucleotideStore& rhs) const
  {
    return d_storage < rhs.d_storage; // XXX ONLY WORKS FOR EQUAL LENGTH, NO CURVAL
  }
  
  static char getVal(char c);
private:
  int bitpos{0};
  std::string d_storage;
  uint8_t d_curval{0};
};

std::ostream& operator<<(std::ostream& os, const NucleotideStore& ns);
std::ostream& operator<<(std::ostream& os, const NucleotideStore::Delta& delta);
