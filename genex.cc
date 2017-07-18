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
  const unsigned int d_hashsize=1<<24;

  uint32_t count(const NucleotideStore& ns, const ReferenceGenome& rg) const;
  vector<pair<uint32_t,bool>> getPositions(const NucleotideStore& ns, const ReferenceGenome& rg, uint32_t before=std::numeric_limits<uint32_t>::max()) const;
  void add(const NucleotideStore& ns, uint32_t pos);
  
} g_hashes;

NucleotideStore g_allA, g_allC, g_allG, g_allT;

void HashCollector::add(const NucleotideStore& stretch, uint32_t pos)
{
  //  
  //  return;
  // CG AT
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

void doKMERCount(const ReferenceGenome& rg)
{
  //  uint32_t uniques=0, reposWin=0;
  std::map<uint32_t, uint32_t> popcounts;
  cout<<"Starting k-mer count"<<endl;

  struct Count
  {
    uint32_t canon{0};
    uint32_t anticanon{0};
    uint32_t count() const { return canon + anticanon; }
  };
  std::vector<Count> lcounts;
  lcounts.resize(std::numeric_limits<uint32_t>::max());

  for(uint32_t pos = 0 ; pos < rg.numNucleotides(); ++pos) {
    if(!(pos % 1024000)) {
      cout<<"\rNow at "<< (100.0*pos / rg.numNucleotides())<<"%";
      cout.flush();
    }
    try {
      auto stretch=rg.getRange(pos, g_unitsize);
      if(!stretch.isCanonical()) {
	stretch = stretch.getRC();
	lcounts[*(uint32_t*)(stretch.getString().c_str())].anticanon++;
      }
      else
	lcounts[*(uint32_t*)(stretch.getString().c_str())].canon++;
    }
    catch(std::exception& e) {
      cerr<<e.what()<<endl;
    }
  }
  cout<<"\nHave "<<lcounts.size()<<" different kmers"<<endl;

  vector<pair<uint32_t,Count> > top;
  top.reserve(count_if(lcounts.begin(), lcounts.end(), [](const auto& a) {
	return a.canon + a.anticanon > 1;
      }));
  cout<<"Of which "<<top.capacity()<<" had >1 use"<<endl;
  unsigned int p=0;
  for(auto iter = lcounts.begin(); iter != lcounts.end() ; ++iter, ++p) {
    if(iter->count() > 1)
      top.push_back({p, *iter});
  }

  cout<<"Put in vector, nth sorting now"<<endl;
  int lim=65536*16;
  
  nth_element(top.begin(), top.begin() + lim, top.end(), [](const auto& a, const auto& b) {
      return a.second.count() > b.second.count();
    });

  cout<<"Sorting first "<<lim<<" elements"<<endl;
  sort(top.begin(), top.begin() + lim, [](const auto& a, const auto& b) {
      return a.second.count() > b.second.count();
    });

  ofstream topfile("top");
  uint32_t cumul=0;
  for(auto iter = top.begin() ; iter != top.begin() + lim ; ++iter) {
    NucleotideStore ns;
    ns.setString(std::string((const char*)&iter->first, (const char*)&iter->first +4));
    topfile<<ns<<"\t"<<ns.getRC()<<"\t"<<iter->second.canon<<"\t"<<iter->second.anticanon<<"\t"<<cumul<<"\n";
    cumul+=iter->second.count();
  }
#if 0
  
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
  ReferenceGenome rg(argv[1], indexChr);

  cout<<"Done reading genome, have "<<rg.numChromosomes()<<" chromosomes, "<<
    rg.numNucleotides()<<" nucleotides"<<endl;

  //  doKMERCount(rg);

  
  string shortname=argc > 2 ? argv[2] : "CM000673.2"; // "CM000663.2";
  auto chromo=rg.getChromosome(shortname);

  uint32_t emitted=0;

  struct Choice
  {
    uint16_t prenucs{0};   
    uint32_t pos{0};
    uint16_t matchnucs{0};
    bool reverse{false};
    int deltalen{0};      // in BYTES
    uint16_t cost() const // cost in BYTES
    {
      return prenucs/4+sizeof(pos)+sizeof(matchnucs)+deltalen+3;
    }
  };


  for(unsigned int beg=0; beg < chromo->chromosome.size(); ) {
    for(int pos=0; pos < 96; pos += g_unitsize) {
      cout<<chromo->chromosome.getRange(beg+pos, g_unitsize)<<"   ";
    }
    cout<<beg+chromo->offset<<" "<<beg<<"@"<<shortname<<endl;
    //    vector<vector<std::tuple<uint32_t,bool,uint16_t,uint16_t>>> positions;
    vector<vector<Choice>> positions;
    vector<std::future<vector<Choice>>> futures;
    for(int pos=0; pos < 96; pos += g_unitsize) {
      futures.emplace_back(std::async(std::launch::async, [chromo,pos,&rg,beg]() {
	    // the starter
	    auto str = chromo->chromosome.getRange(beg+pos, g_unitsize);
	    // where we can find such things in the first big+pos+chromo->offset bytes
	    auto matches=g_hashes.getPositions(str, rg, beg+pos+chromo->offset);
	    {
	      ostringstream ss;
	      ss<<(char)0x1b<<'['<<1+(pos/g_unitsize)*19<<'G'<<'#'<<matches.size();
	      write(1, ss.str().c_str(), ss.str().size());
	    }
	    vector<Choice> ann;

	    uint16_t dsize=0, bestlen=16;
	    for(const auto& m: matches) {
	      unsigned int t;
	      // try longer and longer bits
	      for(t = 16; t < 8193; t < 128 ? ++t : t*=2) {
		auto longerh = chromo->chromosome.getRange(beg+pos, t);
		NucleotideStore longerm;
		if(!m.second) {
		  longerm = rg.getRange(m.first, t);
		}
		else {
		  if(m.first < t)
		    continue;
		  longerm = rg.getRange(m.first-t+g_unitsize, t).getRC();
		}
		if(longerh == longerm) {
		  dsize=0;
		  bestlen=t;
		  continue;
		}
		auto delta = longerm.getDelta(longerh);
		
		if(delta.size() > t/16) // too many errors
		  break;
		bestlen=t;
		dsize = delta.size();
		if(!dsize && bestlen==8192)
		  break;
	      }
	      ann.push_back({0, m.first, bestlen, m.second, dsize});
	      if(!dsize && bestlen==8192)
		break;
	    }
	    
	    sort(ann.begin(), ann.end(),
		 [](const auto& a, const auto& b) {
		   return std::make_tuple(-a.matchnucs, a.deltalen, a.pos) <
		     std::make_tuple(-b.matchnucs, b.deltalen, b.pos);
		   
		 });
	    {
	      ostringstream ss;
	      ss<<(char)0x1b<<'['<<1+(pos/g_unitsize)*19<<'G'<<'!';
	      write(1, ss.str().c_str(), ss.str().size());
	    }

	    return ann;
	  }));
    }

    for(auto& f : futures) {
      auto ann=f.get();
      	    //	    printf("#%-18lu", matches.size());
	    //fflush(stdout);
      for(auto &a : ann)
	a.prenucs=positions.size()*g_unitsize;

      positions.push_back(ann);
    }
    
    printf("\n");

    for(unsigned int n=0; n < 40; ++n) {
      bool some=false;
      for(auto& v : positions) {
	if(n < v.size()) {
	  char o[20];                            // R                           pos
	  snprintf(o,sizeof(o), "%c%u+%u/%u" , v[n].reverse ? 'R':' ', v[n].pos,
		   // len                deltasize
		 v[n].matchnucs, v[n].deltalen);
	  printf("%-19s", o);
	  some=true;
	}
	else
	  printf("                   ");
      }
      printf("\n");
      if(!some)
	break;
    }
    cout<<endl;

    vector<Choice> choices;
    for(auto& v : positions) {
      if(v.size())
	choices.push_back(v[0]);
    }
    sort(choices.begin(), choices.end(), [](const auto& a, const auto& b) {
	return a.prenucs/4 + a.matchnucs/4-a.cost() >  b.prenucs/4 + b.matchnucs/4-b.cost(); // who gets us furthest?
      });

    
    if(!choices.empty()) {
      auto& pick = choices[0];
      cout<<"Best: pos="<<(pick.reverse ? 'R': ' ')<<pick.pos<<", prenucs="<<pick.prenucs<<", matchnucs="<<pick.matchnucs<<", deltalen="<<pick.deltalen<<", cost="<<pick.cost()<<endl;
      for(const auto& c : choices) {
	cout<<"\tOption: pos="<<(c.reverse ? 'R': ' ')<<c.pos<<", prenucs="<<c.prenucs<<", matchnucs="<<c.matchnucs<<", deltalen="<<c.deltalen<<", cost="<<c.cost()<<endl;
      }

      if(pick.prenucs + pick.matchnucs < pick.cost()*4) { // not worth it
	cout<<"Doing plain"<<endl;
	beg+=96;
	emitted+=96/4;
      }
      else {
	beg+=pick.matchnucs+pick.prenucs;
	emitted+=pick.cost();
      }
    }
    else {
      beg+=96;
      emitted+=96/4;
    }
    
    cout<<"Ratio: "<<100.0*emitted/(beg/4)<<"%"<<endl;
  }
  


  
}

