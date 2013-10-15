#include <stdint.h>
#include <stdio.h>
#include <stdexcept>
#include "geneannotated.hh"
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include "misc.hh"
using namespace std;

GeneAnnotationReader::GeneAnnotationReader(const std::string& fname)
{
  FILE* fp=fopen(fname.c_str(), "r");
  if(!fp)
    throw runtime_error("Unable to open "+fname+" for gene annotation reading");

  string line;
  GeneAnnotation ga;
  while(stringfgets(fp, &line)) {
    if(line[0]=='#') {
      cerr<<"Annotations from: "<<line;
      continue;
    }
    const char* p=strtok((char*)line.c_str(), ",\"");
    int field=0;
    do {
      switch(field) {
      case 0:
	ga.tag=p;
	break;
      case 2:
	ga.startPos=atoi(p);
	break;
      case 3:
	ga.stopPos=atoi(p);
	break;
      case 4:
	ga.strand = (*p=='+');
      case 5:
	ga.name=p;
	break;
      }
      field++;
    }while((p=strtok(0, ",\"")));
    d_gas.push_back(ga);
  }
  sort(d_gas.begin(), d_gas.end());
}

vector<GeneAnnotation> GeneAnnotationReader::lookup(uint64_t pos)
{
  vector<GeneAnnotation> ret;
  GeneAnnotation ga;
  ga.startPos = pos;
  gas_t::const_iterator iter =  lower_bound(d_gas.begin(), d_gas.end(), ga);
  //  cout<<"pos: "<<pos<<endl;
  if(iter != d_gas.end()) {

    while(iter != d_gas.begin() && pos < iter->stopPos) {
      //      cout<<iter->name<<" going back!!"<<endl;
      --iter;
    }
    iter++;
    for(; iter != d_gas.end() && iter->startPos < pos && iter->stopPos > pos; ++iter) {
      //      cout << iter->startPos <<", "<<iter->stopPos<<": "<<iter->name<<endl;
      ret.push_back(*iter);
    }
  }
  return ret;
}
