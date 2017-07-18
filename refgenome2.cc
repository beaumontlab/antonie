#include "refgenome2.hh"
#include <stdexcept>
#include <thread>
#include <iostream>
#include "misc.hh"

using namespace std;

NucleotideStore ReferenceGenome::getRange(uint32_t offset, uint32_t len) const
{
  auto iter=std::upper_bound(d_lookup.begin(), d_lookup.end(), offset, [](uint32_t offset, const auto& b) {
      return offset< b->offset;
    });
  
  if(iter == d_lookup.end())
    throw std::range_error("Could not find chromosome for offset "+std::to_string(offset)+" and length "+std::to_string(len));
  --iter;
  if((*iter)->offset <= offset && offset < (*iter)->offset + (*iter)->chromosome.size())
    return (*iter)->chromosome.getRange(offset - (*iter)->offset, len);
  else
    throw std::range_error("Could not find chromosome for offset "+std::to_string(offset)+" and length "+std::to_string(len));
}

ReferenceGenome::ReferenceGenome(const boost::string_ref& fname, std::function<void(ReferenceGenome::Chromosome*, std::string)> idx) : d_fname(fname)
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
	running.emplace_back(idx, chromosome, name);
      }

      string fullname=line+1;
	
      char* spacepos=strchr(line+1, ' ');
    
      if(spacepos)
	*spacepos=0;
      name=line+1;

      if(chromosome)
	seenSoFar += chromosome->chromosome.size();
      d_genome[name].offset = seenSoFar;
      d_genome[name].fullname = fullname;
            
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
    running.emplace_back(idx, chromosome, name);
  }

  fclose(fp);

  for(const auto& c : d_genome) {
    d_lookup.push_back(&c.second);
  }
  sort(d_lookup.begin(), d_lookup.end(), [](const auto& a, const auto& b) {
      return a->offset < b->offset;
    });
       

  
  cout<<"Done reading, awaiting threads"<<endl;
  for(auto& r: running) 
    r.join();
}
