#include "refgenome2.hh"
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
#include <boost/dynamic_bitset.hpp>
#include <unordered_map>

using namespace std;

struct Node
{
  std::string parent;
  vector<std::string> children;
  std::string type;
  std::string tag;
};

std::unordered_map<string, Node> nodes;

int main(int argc, char **argv)
{
  if(argc < 3) {
    cerr<<"Syntax: gfflookup annotations.gff refgenome.fna offset1 [offset2]"<<endl;
    return EXIT_FAILURE;
  }

  GeneAnnotationReader gar(argv[1]);
  cout<<"Done with annotations, got "<<gar.size()<<" of them"<<endl;

  unsigned int bytes=0;
  for(const auto& chromo : {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14",
        "15", "16", "17", "18", "19", "20", "21", "X", "Y", "MT"}) {
    auto anns = gar.getAll(chromo);
    //    map<std::string, Gene> genes;
    boost::dynamic_bitset cdsset;
    for(const auto& a : anns) {
      if(a.type=="chromosome") {
        cdsset.resize(a.stopPos+1);
        break;
      }
    }
    int geneCount{0};
    for(const auto& a : anns) {
      if(a.type == "CDS" || a.type=="exon") { //a.type.get().find("RNA") != string::npos ) {
        for(auto i = a.startPos; i != a.stopPos; ++i)
          cdsset.set(i, true);
      }
      if(a.type=="gene")
        geneCount++;
    }
    cout<<"Total CDS/*RNA length for chromosome "<<chromo<<": "<<cdsset.count()<< " out of "<<cdsset.size()<<", " << (100.0*cdsset.count()/cdsset.size())<<"%, "<<geneCount<<" genes"<<endl;
    bytes += cdsset.count()/4;
  }
  cout<<"Total bytes: "<<bytes<<endl;
  return 0;
  ofstream dot("dot");
  dot<<"digraph D {"<<endl;
  for(const auto& r : gar.lookup("2",162142882, 162152404 )) {
    cout <<"\t"<<r.type<<" "<<r.tag<<" "<<(r.startPos-162142882)<<" - " <<(r.stopPos-162142882) <<", ID="<<r.id<<", Parent="<<r.parent<<", strand: "<<(r.strand ? "+" : "-")<<endl;
    string id;
    if(!r.id.get().empty())
      id = r.id;
    else if(!r.tag.empty())
      id = r.type.get() +" " + r.tag;

    if(!r.parent.get().empty() && !id.empty())
      dot << "\""<< id <<"\" -> \"" << r.parent <<"\""<<endl;
    //    if(r.type=="CDS") 
    //      cout<<"\t\t"<<rg.getChromosome("11")->chromosome.getRange(r.startPos, r.stopPos-r.startPos).getRC()<<endl;
  }

  dot <<"}"<<endl;
  return 0;
  //  for(const auto& chromosome : gar.getChromosomes())
  //    cout<<"\t'"<<chromosome<<"'\n";
  ReferenceGenome rg(argv[2]);
  for(int n =0 ; n < 3; ++n) {
    int p = random() % 1000000;
    auto res = gar.lookup("11", p);
    cout<<"For position "<<p<<", found: "<<endl;
    for(const auto& r: res)
      cout <<"\t"<<r.type<<" "<<r.tag<<" "<<r.startPos<<" - " <<r.stopPos <<", ID="<<r.id<<", Parent="<<r.parent<<", strand: "<<(r.strand ? "+" : "-")<<endl;
  }
  

#if 0
  for(int n = 3; n < argc; ++n) {
    for(const auto& ga : gar.lookup("1", atoi(argv[n]))) {
      cout<<atoi(argv[n])<<'\t'<<ga.startPos<<" - "<<ga.stopPos<<'\t'<<ga.name<<'\t'<<ga.tag<< '\t'<<(ga.strand ? '+' : '-')<<'\t'<<(ga.gene ? "gene" : "") << endl;
      continue;


      if(ga.type=="gene") {
      	string gene = rg.snippet(ga.startPos, ga.stopPos+1);
      	if(!ga.strand)
      		reverseNucleotides(&gene);
	cout<<gene<<endl;
      	for(unsigned int n=0; n < gene.size() ; n+=3) {
      		cout<<DNAToAminoAcid(gene.c_str()+n);
      	}
      	cout<<endl;
      }
    }

    
  }
#endif
}
