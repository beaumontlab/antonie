#include "refgenome2.hh"

#include <iostream>
#include "misc.hh"
#include <set>
#include <algorithm>
#include "dnamisc.hh"
#include <boost/lexical_cast.hpp>
#include <thread>
#include <fstream>
#include "nucstore.hh"
#include <atomic>
#include <mutex>
#include <fstream>
#include <boost/container/small_vector.hpp>
#include <future>
#include <sstream>
#include <sys/prctl.h>

using namespace std;

/* Hash all segments in their canonical form 
   Store a vector where these hashes appear
   To make the uniqueness stats, go through all those positions to check the _actual_ different k-mers they represent
   Emit their individual counts

   Memory usage: we end up storing a 4 byte position for each of 3.5 billion k-mers
*/



struct HashStat
{
  std::mutex* m;
  boost::container::small_vector<uint32_t,2> pos;
};

unsigned int g_unitsize=16;

class HashCollector
{
public:
  HashCollector()
  {
    cout<<"Creating the hashbins"<<endl;
    d_hashes.resize(d_hashsize);
    for(auto & h : d_hashes)
      h.m = new std::mutex;
    cout<<"Done"<<endl;
  }

  ~HashCollector()
  {
    for(auto & h : d_hashes)
      delete h.m;
  }
  vector<HashStat> d_hashes;
  const unsigned int d_hashsize=1<<27;

  uint32_t count(const NucleotideStore& ns, const ReferenceGenome& rg) const;
  vector<pair<uint32_t,bool>> getPositions(const NucleotideStore& ns, const ReferenceGenome& rg, uint32_t before=std::numeric_limits<uint32_t>::max()) const;
  void add(const NucleotideStore& ns, uint32_t pos);
  
} g_hashes;

NucleotideStore g_allA, g_allC, g_allG, g_allT;

void HashCollector::add(const NucleotideStore& stretch, uint32_t pos)
{
  uint32_t h;
  if(!stretch.isCanonical())
    h = stretch.getRC().hash() % d_hashsize;
  else
    h = stretch.hash() % d_hashsize;
    
  //  cout<<"Storing '"<<stretch<<"' at pos, h="<<h<<endl;
  std::lock_guard<std::mutex> l(*d_hashes[h].m);
  if(stretch==g_allA || stretch==g_allC || stretch == g_allG || stretch == g_allT)
    if(d_hashes[h].pos.size() >= 100)
      return;
  d_hashes[h].pos.push_back(pos);  
}

uint32_t HashCollector::count(const NucleotideStore& stretch, const ReferenceGenome& rg) const
{
  uint32_t h;

  if(stretch.isCanonical()) {
    h = stretch.hash() % d_hashsize;
  }
  else
    h = stretch.getRC().hash() % d_hashsize;

  uint32_t ret=0;
  std::lock_guard<std::mutex> l(*d_hashes[h].m);
  
  for(const auto& e : d_hashes[h].pos) {
    auto cmp = rg.getRange(e, g_unitsize);
    if(cmp==stretch)
      ++ret;
    else if(cmp.getRC() == stretch) {
      ++ret;
    }
  }
  return ret;
}

vector<pair<uint32_t,bool>> HashCollector::getPositions(const NucleotideStore& stretch, const ReferenceGenome& rg, uint32_t before) const
{
  vector<pair<uint32_t,bool>> ret;

  uint32_t h;

  if(stretch.isCanonical()) {
    h = stretch.hash() % d_hashsize;
  }
  else
    h = stretch.getRC().hash() % d_hashsize;


  std::lock_guard<std::mutex> l(*d_hashes[h].m);
  //  cout<<"Lookup "<<stretch<<", h="<<h<<", have "<<d_hashes[h].pos.size()<<" candidates"<<endl;
  for(const auto& e : d_hashes[h].pos) {
    if(e >= before)
      continue;
    auto cmp = rg.getRange(e, g_unitsize);
    if(cmp==stretch) 
      ret.push_back({e,false});
    else if (cmp.getRC()==stretch)
      ret.push_back({e,true});
    else {
      //      cout<<"Candidate "<<cmp<<" matched neither way"<<endl;
    }
  }
  return ret;
}

void indexChr(ReferenceGenome::Chromosome* chromosome, std::string name)
{
  if(chromosome->fullname.find("primary")==string::npos) {
    cout<<"Not indexing "<<chromosome->fullname<<endl;
    return;
  }
  
  prctl(PR_SET_NAME, string("Indexing "+name).c_str());
  auto size=chromosome->chromosome.size();
  cout<<"Starting index of '"<<name<<"' with "<<size<<" nucleotides"<<endl;
  for(size_t pos =0; pos < size - g_unitsize; ++pos) {
    auto stretch = chromosome->chromosome.getRange(pos,  g_unitsize);
    g_hashes.add(stretch, chromosome->offset+pos);

  }
  cout<<"Done with index of '"<<name<<"' with "<<size<<" nucleotides"<<endl;
}


namespace std {
    template <>
    struct hash<NucleotideStore> {
        size_t operator () (const NucleotideStore& ns) const { return ns.hash(); }
    };
}


template<typename T, size_t X, size_t Y>
struct Matrix
{
  Matrix()
  {
    values = new T[X*Y];
    for(size_t x=0 ; x < X; ++x)
      for(size_t y=0 ; y < Y; ++y)
	values[x + X*y]=0;
  }
  ~Matrix()
  {
    delete[] values;
  }
  T* values; 
  T& operator()(size_t x, size_t y) {
    return values[x +X*y];
  }

  size_t maxX() const { return X-1; } 
  size_t maxY() const { return Y-1; } 
  
};

Matrix<atomic<uint32_t>, 4000, 4000> g_m;

void doKMERMap(const ReferenceGenome& rg)
{
  auto numnucs= rg.numNucleotides();

  atomic<uint32_t> sofar{0};
  auto f=[&sofar,&rg,&numnucs]() {
    for(uint32_t pos = sofar++ ; pos < numnucs; pos = sofar++) {
      unsigned int xpos =  g_m.maxX() * pos / numnucs;
      unsigned int ypos;

      try {
	auto stretch=rg.getRange(pos, g_unitsize);
	auto matches=g_hashes.getPositions(stretch, rg);
	for(const auto& m : matches) {
	  ypos = g_m.maxY() * m.first/numnucs;
	  g_m(xpos, ypos)++;
	}
      }
      catch(std::exception& e) {
	cerr<<e.what()<<endl;
      }
    }
  };

  vector<std::thread> running;
  for(int n=0; n < 16; ++n)
    running.emplace_back(f);

  while(sofar < numnucs) {
    cout<<sofar<<endl;
    sleep(30);
    ofstream plot("plot");
    for(size_t x=0; x <= g_m.maxX(); ++x)
      for(size_t y=0; y <= g_m.maxY(); ++y)
	plot<<x<<"\t"<<y<<"\t"<<g_m(x,y)<<"\n";
  }
  
  for(auto& r : running)
    r.join();

}


int main(int argc, char**argv)
{
  for(unsigned int n=0; n < g_unitsize;++n) {
    g_allA.append('A');
    g_allC.append('C');
    g_allG.append('G');
    g_allT.append('T');
  }
  
  cout<<"Start reading genome"<<endl;
    
  if(argc < 2) {
    cerr<<"Syntax: genex reference.fasta"<<endl;
    return EXIT_FAILURE;
  }
  ReferenceGenome rg(argv[1], indexChr);

  cout<<"Done reading genome, have "<<rg.numChromosomes()<<" chromosomes, "<<
    rg.numNucleotides()<<" nucleotides"<<endl;

  doKMERMap(rg);
  
}

