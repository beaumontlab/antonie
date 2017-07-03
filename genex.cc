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
#include <boost/container/static_vector.hpp>

using namespace std;

class ReferenceGenome
{
public:
  ReferenceGenome(const boost::string_ref& fname);

private:
  string d_fname;
  map<string,NucleotideStore> d_chromosomes;
};


std::atomic<uint32_t>* g_counts;

std::mutex countlock;
void indexChr(NucleotideStore* chromosome, std::string name)
{
  std::unordered_map<uint32_t, uint32_t> counts;
  auto size=chromosome->size();
  cout<<"Starting index of '"<<name<<"' with "<<size<<" nucleotides"<<endl;
  for(size_t pos =0; pos < size; ++pos) {
    auto stretch = chromosome->getRange(pos,  16);
    g_counts[stretch.hash()]++;
  }
}

ReferenceGenome::ReferenceGenome(const boost::string_ref& fname) : d_fname(fname)
{
  FILE* fp = fopen(d_fname.c_str(), "rb");
  if(!fp)
    throw runtime_error("Unable to open reference genome file '"+d_fname+"'");

  
  char line[256]="";
  string name;
  NucleotideStore* chromosome=0;

  vector<std::thread> running;
  while(fgets(line, sizeof(line), fp)) {
    chomp(line);

    if(line[0] == '>') {
      if(chromosome) {
	running.emplace_back(indexChr,chromosome, name);
      }

      string fullname=line+1;
      
      char* spacepos=strchr(line+1, ' ');
    
      if(spacepos)
	*spacepos=0;
      name=line+1;

      chromosome = &d_chromosomes[name];

      cout<<"Reading chromosome "<<name<<endl;
    }
    else if(chromosome) {
      try {
	chromosome->append(line);
      }
      catch(std::exception& e) {
	cerr<<"Problem storing line "<<line<<endl;
      }
    }
	  
  }
  fclose(fp);
  cout<<"Done reading, awaiting threads"<<endl;
  for(auto& r: running)
    r.join();
}


// stitcher fasta startpos fastq fastq
int main(int argc, char**argv)
{
  g_counts = new std::atomic<uint32_t>[std::numeric_limits<uint32_t>::max()];
  /*
  for(uint64_t pos=0; pos < std::numeric_limits<uint32_t>::max();++pos) {
    g_counts[pos]=0;
  }
  */
  cout<<g_counts[100]<<endl;
  cout<<"Start reading genome"<<endl;
    
  if(argc < 2) {
    cerr<<"Syntax: genex reference.fasta"<<endl;
    return EXIT_FAILURE;
  }
  ReferenceGenome rg(argv[1]);

  cout<<"Done reading genome"<<endl;
  
  uint32_t filled=0, total=0;
  vector<pair<uint32_t, uint32_t>> matches;

  for(uint64_t pos=0; pos < std::numeric_limits<uint32_t>::max();++pos) {
    if(g_counts[pos])
      ++filled;
    total+= g_counts[pos];
    matches.push_back({pos, g_counts[pos]});
  }

  std::sort(matches.begin(), matches.end(), [](const auto&a , const auto& b) {
      return a.second > b.second;
    });
  
  ofstream hashes("hashes");
  for(auto iter = matches.rbegin(); iter != matches.rend(); ++iter) {
    hashes<<iter->first<<"\t"<<iter->second<<"\n";
  }

  
}

