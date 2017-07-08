#include "refgenome.hh"
#include "geneannotated.hh"
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
#include <sys/prctl.h>

using namespace std;

/* Hash all segments in their canonical form 
   Store a vector where these hashes appear
   To make the uniqueness stats, go through all those positions to check the _actual_ different k-mers they represent
   Emit their individual counts

   Memory usage: we end up storing a 4 byte position for each of 3.5 billion k-mers
*/

class ReferenceGenome
{
public:
  ReferenceGenome(const boost::string_ref& fname);

  string d_fname;
  struct Chromosome
  {
    uint32_t offset;
    NucleotideStore chromosome;
  };
  NucleotideStore getRange(uint32_t offset, uint32_t len) const;
  const Chromosome* getChromosome(const std::string& name) const
  {
    if(!d_genome.count(name))
      return 0;
    auto str=d_genome.find(name);
    return &str->second;
  }
  uint32_t numChromosomes()
  {
    return d_genome.size();
  }

  uint32_t numNucleotides()
  {
    uint32_t ret=0;
    for(const auto& g : d_genome)
      ret+=g.second.chromosome.size();
    return ret;
  }

private:
  
  map<string,Chromosome> d_genome;
};

NucleotideStore ReferenceGenome::getRange(uint32_t offset, uint32_t len) const
{
  for(const auto& c : d_genome) {
    if(c.second.offset <= offset && offset < c.second.offset + c.second.chromosome.size())
      return c.second.chromosome.getRange(offset - c.second.offset, len);
  }
  throw std::range_error("Could not find chromosome for offset");
}

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
  const unsigned int d_hashsize=1<<28; 

  uint32_t count(const NucleotideStore& ns, const ReferenceGenome& rg) const;
  vector<pair<uint32_t,bool>> getPositions(const NucleotideStore& ns, const ReferenceGenome& rg) const;
  void add(const NucleotideStore& ns, uint32_t pos);
  
} g_hashes;

NucleotideStore g_allA, g_allC, g_allG, g_allT;

void HashCollector::add(const NucleotideStore& stretch, uint32_t pos)
{
  if(stretch==g_allA || stretch==g_allC || stretch == g_allG || stretch == g_allT)
    return;
  // CG AT
  uint32_t h;
  if(stretch.get(0)=='G' || stretch.get(0)=='T') 
    h = stretch.getRC().hash() % d_hashsize;
  else
    h = stretch.hash() % d_hashsize;
    
  //  cout<<"Storing '"<<stretch<<"' at pos, h="<<h<<endl;
  std::lock_guard<std::mutex> l(*d_hashes[h].m);
  d_hashes[h].pos.push_back(pos);  
}

uint32_t HashCollector::count(const NucleotideStore& stretch, const ReferenceGenome& rg) const
{
  uint32_t h;

  if(stretch.get(0)=='G' || stretch.get(0)=='T')  {
    h = stretch.getRC().hash() % d_hashsize;
  }
  else
    h = stretch.hash() % d_hashsize;

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

vector<pair<uint32_t,bool>> HashCollector::getPositions(const NucleotideStore& stretch, const ReferenceGenome& rg) const
{
  vector<pair<uint32_t,bool>> ret;

  uint32_t h;

  if(stretch.get(0)=='G' || stretch.get(0)=='T')  {
    h = stretch.getRC().hash() % d_hashsize;
  }
  else
    h = stretch.hash() % d_hashsize;


  std::lock_guard<std::mutex> l(*d_hashes[h].m);
  //  cout<<"Lookup "<<stretch<<", h="<<h<<", have "<<d_hashes[h].pos.size()<<" candidates"<<endl;
  for(const auto& e : d_hashes[h].pos) {
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
  prctl(PR_SET_NAME, string("Indexing "+name).c_str());
  auto size=chromosome->chromosome.size();
  cout<<"Starting index of '"<<name<<"' with "<<size<<" nucleotides"<<endl;
  for(size_t pos =0; pos < size - g_unitsize; ++pos) {
    auto stretch = chromosome->chromosome.getRange(pos,  g_unitsize);
    g_hashes.add(stretch, chromosome->offset+pos);

  }
  cout<<"Done with index of '"<<name<<"' with "<<size<<" nucleotides"<<endl;
}


ReferenceGenome::ReferenceGenome(const boost::string_ref& fname) : d_fname(fname)
{
  FILE* fp = fopen(d_fname.c_str(), "rb");
  if(!fp)
    throw runtime_error("Unable to open reference genome file '"+d_fname+"'");

  
  char line[256]="";
  string name;
  ReferenceGenome::Chromosome* chromosome=0;

  vector<std::thread> running;
  uint32_t seenSoFar=0;
  while(fgets(line, sizeof(line), fp)) {
    chomp(line);

    if(line[0] == '>') {
      if(chromosome) {
	running.emplace_back(indexChr, chromosome, name);
      }

      string fullname=line+1;
      
      char* spacepos=strchr(line+1, ' ');
    
      if(spacepos)
	*spacepos=0;
      name=line+1;

      if(chromosome)
	seenSoFar += chromosome->chromosome.size();
      d_genome[name].offset = seenSoFar;
            
      chromosome = &d_genome[name];

      cout<<"Reading chromosome "<<name<<endl;
    }
    else if(chromosome) {
      try {
	chromosome->chromosome.append(line);
      }
      catch(std::exception& e) {
	cerr<<"Problem storing line "<<line<<endl;
      }
    }
	  
  }
  if(chromosome) {
    running.emplace_back(indexChr, chromosome, name);
  }

  fclose(fp);
  cout<<"Done reading, awaiting threads"<<endl;
  for(auto& r: running) 
    r.join();
}

/* 
   MOSTLY WRONG:

   Idea: make dictionary of all K-mers, ordered from most used to least used.
   K-mers btw are only in canonical order (reverse complimented so they start with 'C' or 'A'.

   For k=16 (4 bytes),  this gives us a 500M unique K-mers perhaps. 
   If we allow an edit distance of 1 between K-mers, it is likely we can reduce this to 200M unique k-mers.




   To express any part of the genome, we need to know which of the 200M we want to use, if we need to reverse it,
   and if we need to edit it somehow. This requires a few more bits of information.

   For each k-mer we store where it can first be found in the full.

   To express the genome, we have two modes: actual nucleotides or references.
   We shift from one mode to the other. References are only allowed to be backwards and consist of:
      Location (using variable length encoding) - 30 bits likely
      Reverse complemented or not (1 bit)
      Edit instruction (4 bits where the modification is, 2 bits for new codon)

   To decompress, allocate the whole genome in memory and build it up

 */


// stitcher fasta startpos fastq fastq
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
  ReferenceGenome rg(argv[1]);

  cout<<"Done reading genome, have "<<rg.numChromosomes()<<" chromosomes, "<<
    rg.numNucleotides()<<" nucleotides"<<endl;

  auto chromo=rg.getChromosome("CM000665.2" /*"CM000663.2"*/);

  for(unsigned int beg=0; beg < chromo->chromosome.size(); beg += 96) {
    cout<<chromo->chromosome.getRange(beg, 96)<<" "<<beg+chromo->offset<<" "<<beg<<"@chromosome"<<endl;
    vector<vector<std::tuple<uint32_t,bool,uint16_t>>> positions;
    for(int pos=0; pos < 96; pos += g_unitsize) {
      auto str = chromo->chromosome.getRange(beg+pos, g_unitsize);
      auto matches=g_hashes.getPositions(str, rg);


      vector<std::tuple<uint32_t,bool,uint16_t>> ann;
      auto longerh = chromo->chromosome.getRange(beg+pos, 256);
      for(const auto& m: matches) {
	auto longerm = rg.getRange(std::get<0>(m), 256);
	ann.push_back(std::make_tuple(std::get<0>(m), std::get<1>(m), longerh.fuzOverlap(longerm, 16)));
      }

      sort(ann.begin(), ann.end(),
	   [](const auto&a, const auto& b) {
	     return std::make_tuple(-std::get<2>(a), std::get<0>(a)) <
	       std::make_tuple(-std::get<2>(b), std::get<0>(b));
	   });

      positions.push_back(ann);
      printf("#%-15lu", positions.rbegin()->size());
    }
    printf("\n");
    for(unsigned int n=0; n < 40; ++n) {
      bool some=false;
      for(auto& v : positions) {
	if(n < v.size()) {
	  printf("%c%-11u+%-3u" , std::get<1>(v[n]) ? 'R':' ', std::get<0>(v[n]),
		 std::get<2>(v[n]));
	  some=true;
	}
	else
	  printf("                ");
      }
      printf("\n");
      if(!some)
	break;
    }
    cout<<endl;
  }
  
  uint32_t uniques=0, reposWin=0;
  std::map<uint32_t, uint32_t> popcounts;
  int hcount=0;

  cout<<"Starting k-mer count"<<endl;

  for(const auto& h : g_hashes.d_hashes) {
    if(!(hcount % 102400)) {
      cout<<"\rNow at "<<hcount*100.0/g_hashes.d_hashsize<<"%";
      cout.flush();
    }
    hcount++;

    map<NucleotideStore, uint32_t> lcounts;

    for(const auto& e: h.pos) {
      const auto stretch=rg.getRange(e, g_unitsize);
	lcounts[stretch]++;

      
      
    }

    for(const auto& l : lcounts) {
      popcounts[l.second]++;
    }
  }
  cout<<"\nHave "<<popcounts.size()<<" different counts"<<endl;
  
  {
    ofstream top("top");
    uint64_t tot=0;
    for(const auto& pc : popcounts) {
      top<<pc.first<<"\t"<<pc.second<<"\n";
      tot+=pc.first*pc.second;
    }
    cout<<"Done emitting popcount, total: "<<tot<<endl;
  }
  
  hcount=0;
  
  for(const auto& h : g_hashes.d_hashes) {
    if(!(hcount % 102400))
      cout<<"Now at "<<hcount*100.0/g_hashes.d_hashsize<<"%\n";
    hcount++;
    map<NucleotideStore, uint32_t> lcounts;
    for(const auto& e: h.pos)
      lcounts[rg.getRange(e, g_unitsize)]++;
    
    for(const auto& l : lcounts) {
      // now have bunch of ranges and their counts
      // if we find one with only one count, see if we can improve it
      if(l.second == 2) {
	uniques++;
	uint32_t bestNewCount=l.second;
	//	uint32_t bestNewHash=0;
	NucleotideStore bestNewStretch;

	auto rev=l.first.getRC();
	if(g_hashes.count(rev, rg)) {
	  reposWin++;
	  cout<<reposWin<<"/"<<uniques<<" found "<<l.first<<" on RC to "<<rev<<endl;
	  continue;
	}
	for(unsigned int n = 0; n < l.first.size(); ++n) {
	  for(char c=0; c < 4; ++c) {
	    NucleotideStore nstretch=l.first;
	    nstretch.set(n, c);
	    auto count = g_hashes.count(nstretch, rg);
	    if(count > bestNewCount) {
	      bestNewCount = count;
	      bestNewStretch = nstretch;
	      goto done;
	    }
	  }
	}
      done:;
	if(bestNewCount > 1) {
	  reposWin++;
	  // we would love to actually look this up based on where we found it
	  cout<<reposWin<<"/"<<uniques<<" "<<l.first<<" change to "<<bestNewStretch<<": "<<bestNewCount<<endl;
	}
      }
    }
  }
  


  cout<<reposWin<<" single hits improved after SNP"<<endl;

  
#if 0
  uint32_t filled=0, total=0;
  vector<pair<uint32_t, uint32_t>> matches;
  uint32_t reposWin=0, uniques=0;

  for(uint64_t pos=0; pos < g_mainsize ;++pos) {
    if(!g_counts[pos].pos.empty()) {
      ++filled;
      total+= g_counts[pos].pos.size();
      matches.push_back({pos, g_counts[pos].pos.size()});
    }
    if(g_counts[pos].pos.size() == 1) {
      if(!(pos%16384)) {
	cout<<reposWin<<" upgrades so far out of "<<pos<<" total of which "<<uniques<<" uniques"<<endl;
      }
    }
  }

  cout<<"Starting sort of "<<matches.size()<<" entries"<<endl;
  std::sort(matches.begin(), matches.end(), [](const auto&a , const auto& b) {
      return a.second > b.second;
    });
  cout<<"Start store of hashes"<<endl;
  ofstream hashes("hashes");
  for(auto iter = matches.begin(); iter != matches.end(); ++iter) {
    hashes<<iter->first<<"\t"<<iter->second<<"\n";
  }
#endif
  
}

