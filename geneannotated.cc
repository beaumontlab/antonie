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
#include <boost/algorithm/string.hpp>
using namespace std;


// 1       ensembl_havana  gene    9234775 9271337 .       +       .       ID=gene:ENSG00000049239;Name=H6PD;biotype=protein_coding;description=hexose-6-phosphate dehydrogenase/glucose 1-dehydrogenase [Source:HGNC Symbol%3BAcc:HGNC:4795];gene_id=ENSG00000049239;logic_name=ensembl_havana_gene;version=12

GeneAnnotationReader::GeneAnnotationReader(const std::string& fname)
{
  if(fname.empty())
    return;

  if(!boost::ends_with(fname, ".gff") && !boost::ends_with(fname, ".gff3")) {
    parseGenBank(fname);
    return;
  }
  FILE* fp=fopen(fname.c_str(), "rb");
  if(!fp)
    throw runtime_error("Unable to open '"+fname+"' for gene annotation reading");

  string line;
  map<std::string, vector<Interval<unsigned int, GeneAnnotation> > > gas;
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
      case 0:
        ga.chromosome = p;
        break;
      case 2:
	ga.type=p;
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
      case 8:
	attributeStr=p;
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
      else if(val.first=="ID")
        ga.id = val.second;
      else if(val.first=="Parent")
        ga.parent = val.second;
    }

    if(ga.type =="gene" || ga.type=="CDS" || ga.type=="cds")
      ga.gene=true;
    if(!ga.tag.empty()) {
      ga.tag = ga.type.get() + ": "+ ga.tag;
    }
    else
      ga.tag=ga.type.get();
    
    gas[ga.chromosome].push_back(Interval<unsigned int, GeneAnnotation>(ga.startPos, ga.stopPos, ga));
  no:;
  }
  for(const auto& ga : gas) {
    vector<Interval<unsigned int, GeneAnnotation>> vec = ga.second;
    IntervalTree<unsigned int, GeneAnnotation> tree(std::move(vec));
   
    d_gas[ga.first]=tree;
  }
}

vector<GeneAnnotation> GeneAnnotationReader::getAll(string_view chromo)
{
  vector<GeneAnnotation> ret;
  d_gas[(string)chromo].visit_all([&ret](const auto& i) {
      ret.push_back(i.value);
    });
  return ret;
}

vector<GeneAnnotation> GeneAnnotationReader::lookup(string_view chromo, uint64_t pos1)
{
  vector<GeneAnnotation> ret;

  auto results = d_gas[(string)chromo].findOverlapping(pos1, pos1);
  for(const auto& res : results) {
    ret.push_back(res.value);
  }
  return ret;
}

vector<GeneAnnotation> GeneAnnotationReader::lookup(string_view chromo, uint64_t pos1, uint64_t pos2)
{
  vector<GeneAnnotation> ret;

  auto results = d_gas[(string)chromo].findContained(pos1, pos2);
  for(const auto& res : results) {
    ret.push_back(res.value);
  }
  return ret;
}


void GeneAnnotationReader::parseGenBank(const std::string& fname)
{
  FILE* fp=fopen(fname.c_str(), "rb");
  if(!fp)
    throw runtime_error("Unable to open '"+fname+"' for gene annotation reading");

  string line;

  //
  while(stringfgets(fp, &line)) {
    if(line.find("FEATURES") == 0)
      break;
  }
  
  string genbank;
  while(stringfgets(fp, &line)) {
    if(!isspace(line[0]))
      break;

    boost::trim_right(line);

    genbank+=line+"\n";
  }
  //  d_gas[""]=parseGenBankString(genbank);
}
