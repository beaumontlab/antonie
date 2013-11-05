#include <stdint.h>
#include <stdio.h>
#include <stdexcept>
#include "geneannotated.hh"
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <map>
#include "misc.hh"
using namespace std;

GeneAnnotationReader::GeneAnnotationReader(const std::string& fname)
{
  if(fname.empty())
    return;

  FILE* fp=fopen(fname.c_str(), "r");
  if(!fp)
    throw runtime_error("Unable to open '"+fname+"' for gene annotation reading");

  string line;
  GeneAnnotation ga;
  while(stringfgets(fp, &line)) {
    if(line[0]=='#') {
      //      cerr<<"Annotations from: "<<line;
      continue;
    }
    const char* p=strtok((char*)line.c_str(), "\t\n");
    int field=0;
    string attributeStr;
    do {
      switch(field) {
      case 8:
	attributeStr=p;
	break;
      case 3:
	ga.startPos=atoi(p);
	break;
      case 4:
	ga.stopPos=atoi(p);
	break;
      case 6:
	ga.strand = (*p=='+');
      case 2:
	ga.type=p;
	break;
      }
      field++;
    } while((p=strtok(0, "\t\n")));
    if(ga.type=="repeat_region")
      continue;

    map<string, string> attributes;
    if((p=strtok((char*)attributeStr.c_str(), ";"))) {
      do {
	const char *e = strchr(p, '=');
	if(e) {
	  attributes[string{p,e}]=e+1;
	}
      }while((p=strtok(0, ";")));
    }
    ga.tag.clear();

    for(const auto& val : attributes) {
      if(val.first=="Note" || val.first=="Name" || val.first=="Product") {
	ga.tag.append(val.second);
	ga.tag.append(" ");
      }
      if(val.first=="genome" && val.second=="chromosome")
	goto no;
    }

    d_gas.push_back(ga);
  no:;
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
