#include "misc.hh"
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
/* chromo  source  type    start   stop    ?       strand  ?       key "val" ; key "val";
   1       havana  gene    11869   14409   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene";
*/

using namespace std;

std::string_view getNextWord(const std::string& line, string::size_type* pos)
{
  
  while(line.at(*pos)==' ' || line.at(*pos)=='\t')
    (*pos)++;
  auto startpos = *pos;
  while(line.at(*pos)!=' ' && line.at(*pos)!='\t')
    (*pos)++;

  return std::string_view(&line[startpos], *pos-startpos);
}


std::string_view getNextQuotedWord(const std::string& line, string::size_type* pos)
{
  
  while(line.at(*pos)==' ' || line.at(*pos)=='\t')
    (*pos)++;

  if(line.at(*pos)!='"')
    throw std::runtime_error("Quoted word was not quoted: " + line.substr(*pos));
  (*pos)++;
  auto startpos = *pos;
  while(line.at(*pos)!='"')
    (*pos)++;

  (*pos)++;
  return std::string_view(&line[startpos], *pos-startpos-1);
}


int main(int arvg, char **argv)
{
  std::string line;
  map<string, map<string, map<int,vector<pair<int,int>>>>> exons;
  map<string, int> genesizes;
  map<string, pair<int, int>> genepos;
  map<string, bool> genestrand;
  ofstream csv("genes.csv");
  ofstream exonscsv("exons.csv");
  
  //csv << gene.first << " "<< gexons.second.size() << " " <<totexonsize <<" " <<genesizes[gene.first]<<"\n";
  csv << "chromo name numexons exonsize size start stop strand\n";
  exonscsv << "chromo name start stop strand\n";
  while(stringfgets(stdin, &line)) {
    if(auto pos = line.find('#'); pos != string::npos)
      line.resize(pos);
    if(line.empty())
      continue;

    string::size_type pos = 0;
    auto chromo=getNextWord(line, &pos);
    auto source=getNextWord(line, &pos);
    auto type=getNextWord(line, &pos);
    auto start = atoi(&getNextWord(line, & pos)[0]);
    auto stop = atoi(&getNextWord(line, & pos)[0]);
    auto unknown1 = getNextWord(line, &pos);
    bool sense = getNextWord(line, &pos)=="+";
    auto unknown2 = getNextWord(line, &pos);

    string_view rest{&line.at(pos), line.size()-pos-1};

    string_view gene_name, gene_biotype;
    int transcript_version=-1;
    while(pos +1  < line.size()) {
      auto key = getNextWord(line, &pos);
      auto val = getNextQuotedWord(line, &pos);
      if(key=="gene_name")
        gene_name = val;
      else if(key=="transcript_version")
        transcript_version = atoi(&val[0]);
      else if(key=="gene_biotype")
        gene_biotype = val;
      pos+=1;
    }
    if(gene_biotype != "protein_coding")
      continue;
                 
    
    if(type=="gene") {
      genesizes[(string)gene_name] = stop - start;
      genepos[(string)gene_name] = make_pair(start, stop);
      genestrand[(string)gene_name] = sense;
    }
    else if(type=="exon") {
      cout<<chromo<<"\t"<<sense<<"\t" <<start<<"\t"<<stop<<"\t"<<rest<<endl;

    // gene_id "ENSG00000223972"; gene_version "5"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"

      exons[(string)chromo][(string)gene_name][transcript_version].push_back({start, stop});
    }
  }
  cout<<"Have "<<exons.size()<<" chromosomes"<<endl;
  for(const auto& chromo : exons) {
    cout<<"Chromosome "<<chromo.first<<" has "<<chromo.second.size()<<" gene names"<<endl;
    for(const auto& gene : chromo.second) {
      cout<<gene.first<<" has "<<gene.second.size()<<" transcripts"<<endl;
      int totexonsize=0;
      if(!gene.second.empty()) {
        auto& gexons = *gene.second.rbegin();
        cout <<"   newest transcript has "<<gexons.second.size()<<" exons (";

        for(const auto& exon : gexons.second) {
          exonscsv << chromo.first << " "<<gene.first <<" "<<exon.first<<" "<<exon.second<<" "<<genestrand[gene.first]<<"\n";
          cout<< (exon.second - exon.first) <<" ";
          totexonsize += exon.second - exon.first;

        }
        cout<<")\n";
        csv << chromo.first<<" "<<gene.first << " "<< gexons.second.size() << " " <<totexonsize <<" " <<genesizes[gene.first]<<" "<<genepos[gene.first].first<<" "<<genepos[gene.first].second<<" "<<genestrand[gene.first]<<"\n";
      }
    }
  }
}
