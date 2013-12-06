#include "refgenome.hh"
#include "geneannotated.hh"
#include <iostream>
using namespace std;

int main(int argc, char**argv)
{
  ReferenceGenome rg(argv[1]);
  //  GeneAnnotationReader gar(argv[2]);

  cout<<rg.snippet(atoi(argv[2]), atoi(argv[3]))<<endl;
}
