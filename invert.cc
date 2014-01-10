#include <iostream>
#include "misc.hh"
using namespace std;

int main(int argc, char**argv)
{
  for(int n = 1 ; n < argc; ++n) {
    string nucs(argv[n]);
    reverseNucleotides(&nucs);
    cout<<nucs<<endl;
  }
}
