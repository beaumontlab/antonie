#include "nucstore.hh"
#include <iostream>

using std::cout;
using std::endl;


/* byte  byte byte
   ACGT  ACGT ACGT */

NucleotideStore NucleotideStore::getRange(size_t pos, size_t len) const
{
  NucleotideStore ret;
  for(; len; pos++, len--)
    ret.append(get(pos));
  return ret;
}

char NucleotideStore::get(size_t pos) const
{
  uint8_t byte;
  if(pos/4 < d_storage.size())
    byte=d_storage.at(pos/4);
  else
    byte=d_curval;

  byte >>= ((pos%4)*2);
  return "ACGT"[byte & 0x3];
}

void NucleotideStore::append(const boost::string_ref& line)
{
  for(const auto& c : line)
    append(c);
}

void NucleotideStore::append(char c)
{
  uint8_t val;
  switch(c) {
  case 0:
  case 'A':
  case 'a':
    val=0;
    break;

  case 1:
  case 'C':
  case 'c':
    val=1;
    break;

  case 2:
  case 'G':
  case 'g':
    val=2;
    break;

  case 3:
  case 'T':
  case 't':
    val=3;
    break;
  case 'N':
    return;
  default:
    throw std::runtime_error("Impossible nucleotide: "+std::string(1, c));
  }
  if(!bitpos) {
    d_curval=0;
  }
  val<<=bitpos;
  d_curval|=val;
  bitpos += 2;
  if(bitpos == 8) {
    d_storage.append(1, d_curval);
    bitpos=0;
  }
    
}

std::ostream& operator<<(std::ostream& os, const NucleotideStore& ns)
{
  for(size_t pos=0; pos < ns.size(); ++pos)
    os<<ns.get(pos);
  return os;
}
