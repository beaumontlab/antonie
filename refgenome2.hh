#pragma once
#include <vector>
#include <map>
#include "nucstore.hh"
#include <functional>

class ReferenceGenome
{
public:
  struct Chromosome
  {
    std::string fullname;
    uint32_t offset;
    NucleotideStore chromosome;
  };

  ReferenceGenome(const boost::string_ref& fname, std::function<void(Chromosome*, std::string)> idx);
 
  std::string d_fname;
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

  uint32_t numNucleotides() const
  {
    if(d_lookup.empty())
      return 0;
    return (*d_lookup.rbegin())->offset + (*d_lookup.rbegin())->chromosome.size();
  }

private:
  
  std::map<std::string,Chromosome> d_genome;
  std::vector<const Chromosome*> d_lookup;
};
