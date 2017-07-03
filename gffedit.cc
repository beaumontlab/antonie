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
  ReferenceChromosome rg(argv[1]);
  FILE* fp =fopen(argv[2], "rb");
  if(!fp)
    throw runtime_error("Unable to open '"+string(argv[2])+"' for reading GFF3: "+string(strerror(errno)));
  
  string line;

  map<string, unsigned int> scounts;
  ofstream newgff(argv[3]);
  dnapos_t startInsert=argc > 4 ? atoi(argv[4]) : 0;
  int shiftInsert=argc > 5 ? atoi(argv[5]) : 0;
  while(stringfgets(fp, &line)) {
    GeneAnnotation ga;
    if(!line.empty() && line[0]=='#') {
      newgff<<line;
      continue;
    }
    const char* p=strtok((char*)line.c_str(), "\t\n");
    int field=0;
    string attributeStr;
    do {
      if(field)
	newgff<<'\t';
      switch(field) {
      case 2:
	ga.type=p;
	newgff<<p;
	break;
      case 3:
	ga.startPos=atoi(p);
	if(ga.startPos > startInsert) 
	  ga.startPos += shiftInsert;
	newgff<<ga.startPos;
	break;
      case 4:
	ga.stopPos=atoi(p);
	if(ga.stopPos > startInsert) 
	  ga.stopPos += shiftInsert;

	newgff<<ga.stopPos;
	break;
      case 6:
	ga.strand = (*p=='+');
	newgff<<p;
	break;
      default:
	newgff<<p;
      }
      field++;

    } while((p=strtok(0, "\t\n")));
    newgff<<'\n';
      
    if(ga.type=="gene") {
      if(ga.strand)
	scounts[rg.snippet(ga.startPos, ga.startPos+3)]++;
      else {
	string codon=rg.snippet(ga.stopPos-2, ga.stopPos+1);
	reverseNucleotides(&codon);
	scounts[codon]++;
      }
    }
  }
  vector<pair<unsigned int, string>> rcounts;
  int totCount=0;
  for(const auto& c : scounts) {
    rcounts.push_back({c.second, c.first});
    totCount+=c.second;
  }
  sort(rcounts.begin(), rcounts.end());
  for(auto c = rcounts.rbegin(); c != rcounts.rend(); ++c) {
    cout<<100.0*c->first/totCount<<"%\t"<<c->second<<'\n';
    if(c-rcounts.rbegin() > 10)
      break;
  }
}
