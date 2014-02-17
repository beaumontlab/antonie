#include "refgenome.hh"
#include "geneannotated.hh"
#include <errno.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include "misc.hh"

using namespace std;

// gffedit fasta gff newgff [insertpos insertlen]
int main(int argc, char **argv)
{
  GeneAnnotationReader gar(argv[1]);
  for(int n = 2; n < argc; ++n) {
    for(const auto& ga : gar.lookup(atoi(argv[n])))
      cout<<atoi(argv[n])<<'\t'<<ga.name<<'\t'<<ga.tag<<endl;
    
  }
}
