#include "16ssearcher.hh"
#include <zlib.h>
#include <stdexcept>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include "fastq.hh"
#include "fastqindex.hh"
#include "stitchalg.hh"
#include <fstream>

extern "C" {
#include "hash.h"

}
#include "misc.hh"
using namespace std;

Search16S::Search16S(const std::string& fname)
{
  d_fp = fopen(fname.c_str(), "r");
  if(!d_fp) 
    throw runtime_error("Unable to open file '"+fname+"' for reading");
}

bool Search16S::get(Entry* entry, uint64_t* offset)
{
  char line[16384];
 
  if(offset)
    *offset=ftell(d_fp);

  if(!fgets(line, sizeof(line),d_fp))
    return false;

  if(line[0]!='>') 
    throw runtime_error("Unable to parse line '"+string(line)+"' as green genes 16s line");

  entry->id = atoi(line+1);
  
  if(!fgets(line, sizeof(line), d_fp))
    return false;
  
  if(line[0]!='C' && line[0]!='A' && line[0]!='G' && line[0]!='T' && line[0]!='K' && line[0]!='N' && line[0]!='M' && line[0]!='B')  // K??
    throw runtime_error("Unable to parse line '"+string(line)+"' as green genes 16s line");

  chomp(line);
  entry->nucs = line;
  return true;
}

void Search16S::seek(uint64_t pos)
{
  fseek(d_fp, pos, SEEK_SET);
}


string nameFromAccessionNumber(const std::string& number)
{
  static map<string, string> s_names;
  if(s_names.count(number))
    return s_names[number];
  
  FILE* fp = popen((string("wget -q -O- 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=docsum&email=bert.hubert@netherlabs.nl&id=")+number+"' | grep Title | cut -f2 -d'>' | cut -f1 -d'<'").c_str(), "r");
  
  if(!fp) 
    throw runtime_error("Unable to open wget pipe: "+string(strerror(errno)));

  string line;
  if(stringfgets(fp, &line)) {
    boost::trim_right(line);
    s_names[number]=line;
  }
  pclose(fp);
  return line;
  
}

map<int, string> GreenGenesToGenbank(const std::string& fname)
{
  FILE* fp = fopen(fname.c_str(), "r");
  if(!fp)
    throw runtime_error("Unable to open "+fname+" for reading: "+string(strerror(errno)));
  
  map<int, string> ret;
  char line[1024];
  int id;
  char* p;
  while(fgets(line, sizeof(line), fp)) {
    if(*line=='#')
      continue;
    id=atoi(line);
    if(!id)
      throw runtime_error("This does not look like a Green Genes accession map!");
    p=strchr(line, '\n');
    if(p) *p =0;
    
    p=strchr(line, '\t');
    if(p) {
      ret[id]=p+1;
    }
  }
  fclose(fp);
  return ret;
};

/* idea - go through 16S database, score for each entry how many hits we find in the FASTQ index
   for the first n index-length chunks. Order the entries on score, and start stitching them & report 
   all >99% matches */

// 16ssearcher gg_13_5.fasta 1.fastq 2.fastq etc 
int main(int argc, char** argv)
{
  if(argc < 3) {
    cerr<<"Syntax: 16ssearcher gg_13_5.fasta gg_13_5_accessions.txt fastq1 [fastq2..]"<<endl;
  }
  auto idmap = GreenGenesToGenbank(argv[2]);
  map<FASTQReader*, unique_ptr<vector<HashedPos> > > fhpos;

  FASTQReader* fqreader;

  for(int f = 3; f < argc; ++f) {
    fqreader = new FASTQReader(argv[f], 33, 0);
    fhpos[fqreader]=indexFASTQ(fqreader, argv[f], 35);
  }

  Search16S s16(argv[1]);

  Search16S::Entry entry;
  uint64_t offset;
  string part;

  vector<pair<int,uint64_t> > scores;
  int scount=0;
  int maxscore=0;

  while(s16.get(&entry, &offset)) {
    vector<int> qscores(entry.nucs.size());
    for(unsigned int n=0; n < entry.nucs.size()-35; ++n) {
      part= entry.nucs.substr(n, 35);
      for(const auto& match : getConsensusMatches(part, fhpos, 35)) {
	if(dnaDiff(match.d_nucleotides, entry.nucs.substr(n)) < 2 ) {
	  for(string::size_type pos = 0 ; pos < match.d_nucleotides.length() && n+pos < qscores.size(); ++pos) {
	    if(match.d_nucleotides[pos] == entry.nucs[n+pos])
	      qscores[n+pos]+=match.d_quality[pos];
	  }
	  break;
        }
      }
      if(n > 20 && !accumulate(qscores.begin()+n-20, qscores.begin()+n, 0))
	break;
    }
    auto sum=accumulate(qscores.begin(), qscores.end(), 0);
    unsigned int qpos;
    for(qpos = 0; qpos < qscores.size(); ++qpos) {
      if(qpos > 20 && !qscores[qpos])
	break;
    }
    if(qpos == qscores.size()) {
      cerr<<"\nMay be a winner, full coverage ("<<sum/entry.nucs.size()<<") on " <<entry.id<<" -> "<<idmap[entry.id]<< ": "<< nameFromAccessionNumber(idmap[entry.id])<<endl;
      scores.push_back({sum/entry.nucs.size(), offset});
      if(sum >= maxscore) {
	maxscore = sum;
	cerr<<" -> Best current guess: "<<entry.id<<" ("<<sum/entry.nucs.size()<<") -> "<<idmap[entry.id]<<endl;
      }
    }
    if(!((++scount)%1000) || sum)
      cerr<<'\r'<<scores.size()<< " potentials out of "<<scount<<" candidates";
  }

  cerr<<"\rHave "<<scores.size()<<" potentials out of "<<scount<<" candidates"<<endl;

  sort(scores.begin(), scores.end());

  for(auto iter = scores.rbegin(); iter != scores.rend(); ++iter) {
    s16.seek(iter->second);
    s16.get(&entry);
    cout<<"Score: "<<iter->first<<", Green Genes ID: "<<entry.id<< " -> "<<idmap[entry.id]<<": "<< nameFromAccessionNumber(idmap[entry.id])<<endl;;
    string found = doStitch(fhpos, entry.nucs.substr(0, 100), entry.nucs.substr(entry.nucs.length()-100), 1.2*entry.nucs.length(), 35, false);
    cout <<"Diff: "<<dnaDiff(found, entry.nucs)<<endl;
  }
}
