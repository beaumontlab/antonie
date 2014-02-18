#include "refgenome.hh"
#include "geneannotated.hh"
#include <errno.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include "dnamisc.hh"
#include <map>
#include <vector>
#include <algorithm>
#include "misc.hh"

using namespace std;

// gffedit fasta gff newgff [insertpos insertlen]
int main(int argc, char **argv)
{
  GeneAnnotationReader gar(argv[1]);
  ReferenceGenome rg(argv[2]);
  for(int n = 3; n < argc; ++n) {
    for(const auto& ga : gar.lookup(atoi(argv[n]))) {
      cout<<atoi(argv[n])<<'\t'<<ga.name<<'\t'<<ga.tag<< '\t'<<(ga.strand ? '+' : '-')<<endl;
      if(ga.type=="gene") {
      	string gene = rg.snippet(ga.startPos, ga.stopPos+1);
      	if(!ga.strand)
      		reverseNucleotides(&gene);
      	cout<<gene<<endl;
      	for(int n=0; n < gene.size() ; n+=3) {
      		cout<<DNAToAminoAcid(gene.c_str()+n);
      	}
      	cout<<endl;
      }
    }

    
  }
}
