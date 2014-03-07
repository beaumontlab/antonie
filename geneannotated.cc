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

  while(stringfgets(fp, &line)) {
    GeneAnnotation ga;
    ga.gene=false;
    if(line[0]=='#') {
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
	break;
      case 2:
	ga.type=p;
	break;
      }
      field++;
    } while((p=strtok(0, "\t\n")));
    //    if(ga.type=="repeat_region")
    //  continue;


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
      if(val.first=="Note" || val.first=="Name" || val.first=="Product" || val.first=="product") {
	ga.tag.append(val.second);
	ga.tag.append(" ");
      }
      if(val.first=="genome" && val.second=="chromosome")
	goto no;
    }

    if(ga.type =="gene" || ga.type=="CDS" || ga.type=="cds")
      ga.gene=true;
    if(!ga.tag.empty()) {
      ga.tag = ga.type+": "+ga.tag;
      d_gas.push_back(ga);
    }

      
  no:;
  }
}

vector<GeneAnnotation> GeneAnnotationReader::lookup(uint64_t pos)
{
  vector<GeneAnnotation> ret;

  for(const auto& a : d_gas) {
    if(a.startPos <= pos && pos <= a.stopPos)
      ret.push_back(a);
  }
  return ret;
}
