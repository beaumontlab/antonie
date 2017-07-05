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

NucleotideStore NucleotideStore::getRC() const
{
  NucleotideStore ret;
  char c;
  for(size_t pos=size(); pos; --pos) {
    c=get(pos-1);
    if(c=='C')
      c='G';
    else if(c=='G')
      c='C';
    else if(c=='A')
      c='T';
    else if(c=='T')
      c='A';
    ret.append(c); // C -> G, A->t
  }
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

bool getBit(uint8_t c, uint8_t bit)
{
  uint8_t mask = 1 << bit;
  return c & mask;
}




static void setBit(uint8_t * val, int bit, bool nval)
{
  uint8_t mask = 1 << bit;
  if(nval) {
    *val |= mask;
  }
  else
    *val &= (~mask);
}


void NucleotideStore::set(size_t pos, char c) 
{
  uint8_t newval=getVal(c);
  if(pos/4 < d_storage.size()) {
    uint8_t val=d_storage.at(pos/4);
    setBit(&val, (pos%4)*2, getBit(newval, 0));
    setBit(&val, (pos%4)*2+1, getBit(newval, 1));
    
    d_storage.at(pos/4) = val;
  }
  else {
    setBit(&d_curval, (pos%4)*2, getBit(newval, 0));
    setBit(&d_curval, (pos%4)*2+1, getBit(newval, 1));
  }

}


void NucleotideStore::append(const boost::string_ref& line)
{
  for(const auto& c : line)
    append(c);
}

char NucleotideStore::getVal(char c)
{
  switch(c) {
  case 0:
  case 'A':
  case 'a':
    return 0;

    
  case 1:
  case 'C':
  case 'c':
    return 1;
    
  case 2:
  case 'G':
  case 'g':
    return 2;
    
  case 3:
  case 'T':
  case 't':
    return 3;
  case 'N':
    return 4;
  }
  throw std::runtime_error("Impossible nucleotide: "+std::string(1, c));
}

void NucleotideStore::append(char c)
{
  uint8_t val=getVal(c);
  if(val == 4) // N
    return;
  
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
