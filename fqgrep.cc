#include <iostream>
#include "misc.hh"
#include "fastq.hh"
using namespace std;

int main(int argc, char**argv)
{
  string search(argv[1]);
  string rsearch(search);
  reverseNucleotides(&rsearch);

  for(int n = 2 ; n < argc; ++n) {
    FASTQReader fqreader(argv[n], 33);
    FastQRead fqr;

    while(fqreader.getRead(&fqr)) {
      auto pos = fqr.d_nucleotides.find(search);
      if(pos != string::npos) {
	cout<<fqr.d_nucleotides.substr(pos)<<endl;
      }
      else if((pos=fqr.d_nucleotides.find(rsearch)) != string::npos) {
	fqr.reverse();
	pos = fqr.d_nucleotides.find(search);
	cout<<fqr.d_nucleotides.substr(pos)<<endl;
      }
    }
  }
}
