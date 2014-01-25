#include "16ssearcher.hh"
#include <zlib.h>
#include <stdexcept>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <boost/foreach.hpp>
#include "fastq.hh"
#include "fastqindex.hh"
#include "stitchalg.hh"

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
  auto idmap = GreenGenesToGenbank("./gg_13_5_accessions.txt");
  map<FASTQReader*, unique_ptr<vector<HashedPos> > > fhpos;

  FASTQReader* fqreader;

  for(int f = 2; f < argc; ++f) {
    fqreader = new FASTQReader(argv[f], 33, 0);
    fhpos[fqreader]=indexFASTQ(fqreader, argv[f], 35);
  }

  Search16S s16(argv[1]);

  Search16S::Entry entry;
  uint64_t offset;
  string part;
  int matchcount=0;
  vector<pair<int,uint64_t> > scores;
  int scount=0;
  int maxscore=0;
  while(s16.get(&entry, &offset)) {
    matchcount=0;
    for(unsigned int n=0; n < entry.nucs.size()-35; ++n) {
      part= entry.nucs.substr(n, 35);
      for(auto match : getConsensusMatches(part, fhpos, 35)) {
	if(!dnaDiff(match.d_nucleotides,entry.nucs.substr(n)))
	  matchcount++;
      }
      if(n==20 && !matchcount)
	break;
      if(n==100 && matchcount < 5)
	break;
    }
    if(matchcount) {
      //      cout << entry.id << ' ' << matchcount << endl;
      scores.push_back({matchcount, offset});
      if(matchcount >= maxscore) {
	maxscore = matchcount;
	cerr<<" -> Best current guess: "<<entry.id<<" ("<<maxscore<<") -> "<<idmap[entry.id]<<endl;
      }
    }
    if(!((++scount)%1000) || matchcount)
      cerr<<'\r'<<scores.size()<< " potentials out of "<<scount<<" candidates";
  }

  cerr<<"\rHave "<<scores.size()<<" potentials out of "<<scount<<" candidates"<<endl;
  sort(scores.begin(), scores.end());


  for(auto iter = scores.rbegin(); iter != scores.rend(); ++iter) {

    cout << iter->second <<": "<<iter->first<<endl;
    s16.seek(iter->second);
    s16.get(&entry);
    cout<<"ID: "<<entry.id<< " -> "<<idmap[entry.id]<<endl;
    cerr<<"ID: "<<entry.id<< " -> "<<idmap[entry.id]<<endl;
    string found = doStitch(fhpos, entry.nucs.substr(0, 100), entry.nucs.substr(entry.nucs.length()-100), 10000, 35);
    cout <<"Diff: "<<dnaDiff(found, entry.nucs)<<endl;
  }
}
