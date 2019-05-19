#include "refgenome2.hh"
#include <errno.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include "dnamisc.hh"
#include <map>
#include <vector>
#include <algorithm>
#include "misc.hh"
#include <boost/dynamic_bitset.hpp>
#include <unordered_map>

using namespace std;


int main(int argc, char **argv)
{
  if(argc < 3) {
    cerr<<"Syntax: gfflookup refgenome.fna chromosome offset1 offset2"<<endl;
    return EXIT_FAILURE;
  }

  ReferenceGenome rg(argv[1]);
  auto gene = rg.getChromosome(argv[2])->chromosome.getRange(atoi(argv[3]), atoi(argv[4]) - atoi(argv[3]));

  cout<<gene<<endl;
  cout<<gene.getRC()<<endl;

  for(int o=0; o < 3; ++o) {
    auto str = gene.toASCII();
    for(unsigned int n=o; n < str.size() ; n+=3) {
      cout<<DNAToAminoAcid(str.c_str()+n);
    }
    cout<<endl;
    str = gene.getRC().toASCII();
    for(unsigned int n=o; n < str.size() ; n+=3) {
      cout<<DNAToAminoAcid(str.c_str()+n);
    }
    cout<<endl;
  }
  // reverseNucleotides(&gene);
  /*
  for(unsigned int n=0; n < gene.size() ; n+=3) {
    cout<<DNAToAminoAcid(gene.c_str()+n);
  }
  */

}
