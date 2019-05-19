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

std::vector<NucleotideStore::Delta> NucleotideStore::getDelta(const NucleotideStore& b, double mispen, double gappen, double skwpen) const
{
  const NucleotideStore& a=*this;
  std::vector<NucleotideStore::Delta> ret;
  unsigned int i,j,k;
  double dn,rt,dg;
  std::string::size_type ia = size(), ib = b.size();
  double *cost[ia+1];
  for(unsigned int n=0; n < ia+1; ++n)
    cost[n]=new double[ib+1];
  cost[0][0] = 0.;
  for (i=1;i<=ia;i++) cost[i][0] = cost[i-1][0] + skwpen;
  for (i=1;i<=ib;i++) cost[0][i] = cost[0][i-1] + skwpen;
  for (i=1;i<=ia;i++) for (j=1;j<=ib;j++) {
      dn = cost[i-1][j] + ((j == ib)? skwpen : gappen);
      rt = cost[i][j-1] + ((i == ia)? skwpen : gappen);
      dg = cost[i-1][j-1] + ((get(i-1) == b.get(j-1))? -1. : mispen);
      cost[i][j] = std::min({dn,rt,dg});
    }
  i=ia; j=ib; k=0;
  while (i > 0 || j > 0) {
    dn = rt = dg = 9.99e99;
    if (i>0) dn = cost[i-1][j] + ((j==ib)? skwpen : gappen);
    if (j>0) rt = cost[i][j-1] + ((i==ia)? skwpen : gappen);
    if (i>0 && j>0) dg = cost[i-1][j-1] +
		      ((a[i-1] == b[j-1])? -1. : mispen);
    if (dg <= std::min(dn,rt)) {
      //      aout[k] = ain[i-1];
      //      bout[k] = bin[j-1];
      bool match=(a[i-1] == b[j-1]);
      //      summary[k++] = (match ? '=' : '!');
      if(match)
	; //ret.matches++;
      else 
	ret.push_back({(i-1), b[j-1], Delta::Action::Replace}); // ret.mismatches++;
      
      i--; j--;
    }
    else if (dn < rt) {
      //      aout[k] = ain[i-1];
      //      bout[k] = ' ';
      // summary[k++] = ' ';		
      if(j==ib)
        ret.push_back({i-1, 0, Delta::Action::Delete});
      else
        ret.push_back({i-1, 0, Delta::Action::Delete});
	; // ret.deletes++;
      i--;
    }
    else {
      //      aout[k] = ' ';
      // bout[k] = bin[j-1];
      //  summary[k++] = ' ';		
      if(i==ia)
	; // ret.skews++;
      else
        ret.push_back({i-1, b[j-1], Delta::Action::Insert});
	; // ret.inserts++;
      j--;
    }
  }
  for (i=0;i<k/2;i++) {
    //    swap(aout[i],aout[k-1-i]);
    // swap(bout[i],bout[k-1-i]);
    // swap(summary[i],summary[k-1-i]);
  }
  // aout.resize(k); bout.resize(k); summary.resize(k);
  for(unsigned int n=0; n < ia+1; ++n)
    delete[] cost[n];
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

size_t NucleotideStore::overlap(const NucleotideStore& rhs) const
{
  size_t pos=0;
  for(; pos < size() && pos < rhs.size() && get(pos)==rhs.get(pos); ++pos)
    ;
  return pos;
}

size_t NucleotideStore::fuzOverlap(const NucleotideStore& rhs, int ratio) const
{
  size_t pos=0, mism=0;
  for(; pos < size() && pos < rhs.size() && mism <= pos/ratio; ++pos) {
    if(get(pos)!=rhs.get(pos))
      mism++;
  }
  return pos;
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
  default:
    return 0;
  }
  throw std::runtime_error("Impossible nucleotide: "+std::string(1, c));
}

void NucleotideStore::append(char c)
{
  uint8_t val=getVal(c);
  
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

std::string NucleotideStore::toASCII() const
{
  std::string ret;
  ret.reserve(size());
  for(size_t pos=0; pos < size(); ++pos)
    ret.append(1, get(pos));
  return ret;
}

std::ostream& operator<<(std::ostream& os, const NucleotideStore& ns)
{
  for(size_t pos=0; pos < ns.size(); ++pos)
    os<<ns.get(pos);
  return os;
}

std::ostream& operator<<(std::ostream& os, const NucleotideStore::Delta& d)
{
  os<<d.pos<<" "<<d.o<<" ";
  if(d.a==NucleotideStore::Delta::Action::Replace)
    os<<"R";
  else
    if(d.a==NucleotideStore::Delta::Action::Delete)
      os<<"D";
  else
    if(d.a==NucleotideStore::Delta::Action::Insert)
      os<<"I";
  
  return os;
}
