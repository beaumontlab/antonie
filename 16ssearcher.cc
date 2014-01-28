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
#include <boost/lexical_cast.hpp>
#include <fstream>
#include <set>

extern "C" {
#include "hash.h"

}
#include "misc.hh"
using namespace std;

Search16S::Search16S(const std::string& fname)
{
  d_linereader = LineReader::make(fname);
}

bool Search16S::get(Entry* entry, uint64_t* offset)
{
  char line[16384];
 
  if(offset)
    *offset = d_linereader->getUncPos();

  if(!d_linereader->fgets(line, sizeof(line)))
    return false;

  if(line[0]!='>') 
    throw runtime_error("Unable to parse line '"+string(line)+"' as green genes 16s line");

  entry->id = atoi(line+1);
  
  if(!d_linereader->fgets(line, sizeof(line)))
    return false;
  
  if(line[0]!='C' && line[0]!='A' && line[0]!='G' && line[0]!='T' && line[0]!='K' && line[0]!='N' && line[0]!='M' && line[0]!='B')  // K??
    throw runtime_error("Unable to parse line '"+string(line)+"' as green genes 16s line");

  chomp(line);
  entry->nucs = line;
  return true;
}

void Search16S::seek(uint64_t pos)
{
  d_linereader->seek(pos);
}


string nameFromAccessionNumber(const std::string& number)
{
  static map<string, string> s_names;
  if(s_names.count(number))
    return s_names[number];
  
  FILE* fp = popen((string("wget -q -O- 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=docsum&tool=antonie&id=")+number+"' | grep Title | cut -f2 -d'>' | cut -f1 -d'<'").c_str(), "r");
  
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
  auto fp = LineReader::make(fname);
  
  map<int, string> ret;
  char line[1024];
  int id;
  char* p;
  while(fp->fgets(line, sizeof(line))) {
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
    return EXIT_FAILURE;
  }
  auto idmap = GreenGenesToGenbank(argv[2]);
  map<FASTQReader*, unique_ptr<vector<HashedPos> > > fhpos;

  FASTQReader* fqreader;

  for(int f = 3; f < argc; ++f) {
    fqreader = new FASTQReader(argv[f], 33, 0);
    fhpos[fqreader]=indexFASTQ(fqreader, argv[f], 35);
  }

  Search16S s16(argv[1]);


  int scount=0;
  int maxscore=0;

  struct Candidate
  {
    int score;
    Search16S::Entry entry;
    vector<FastQRead> reads;
  };
  vector<Candidate> candidates;
  Candidate candidate;
  while(s16.get(&candidate.entry)) {
    candidate.reads.clear();
    vector<int> qscores(candidate.entry.nucs.size());
    unsigned int n;

    for( n=0; n < candidate.entry.nucs.size()-35; ++n) {
      string part=candidate.entry.nucs.substr(n);
      for(const auto& match : getConsensusMatches(part, fhpos, 35)) {
	if(dnaDiff(match.d_nucleotides, part) < 2 ) {
	  candidate.reads.push_back(match);
	  for(string::size_type pos = 0 ; pos < match.d_nucleotides.length() && n+pos < qscores.size(); ++pos) {
	    if(match.d_nucleotides[pos] == candidate.entry.nucs[n+pos])
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
      ofstream cov(boost::lexical_cast<string>(candidates.size())+".cov");
      for(qpos = 0; qpos < qscores.size(); ++qpos) {
	cov<<qpos<<'\t'<< qscores[qpos]<<'\n';
      }
      
      cerr<<"\nFull coverage ("<<sum/candidate.entry.nucs.size()<<") on " <<
	candidate.entry.id<<" -> "<<idmap[candidate.entry.id]<< ": "<< nameFromAccessionNumber(idmap[candidate.entry.id])<<endl;
      candidate.score=sum/candidate.entry.nucs.size();

      candidates.push_back(candidate);
	  
      if(sum >= maxscore) {
	maxscore = sum;
	cerr<<" -> Best current guess: "<<candidate.entry.id<<" ("<<sum/candidate.entry.nucs.size()<<") -> "<<idmap[candidate.entry.id]<<endl;
      }
    }
    if(!((++scount)%1000) || sum)
      cerr<<'\r'<<candidates.size()<< " potentials out of "<<scount<<" 16S entries";
  }

  cerr<<"\rHave "<<candidates.size()<<" potentials out of "<<scount<<" 16S entries"<<endl;

  sort(candidates.begin(), candidates.end(), 
       [](const Candidate& a, const Candidate& b) { 
	 return a.score < b.score;
       });
  
  cerr<<"Done sorting"<<endl;
  set<FastQRead> allReads;
  for(const auto& c : candidates) {
    allReads.insert(c.reads.begin(), c.reads.end());
  }
  cout<<"Total reads contributing to 16S matches: "<<allReads.size()<<endl;
  for(auto iter = candidates.rbegin(); iter != candidates.rend(); ++iter) {
    cout<<"Score: "<<iter->score<<" ("<<iter->reads.size()*100.0/allReads.size()<<"% of reads), Green Genes ID: "<<iter->entry.id<< " -> "<<idmap[iter->entry.id]<<": "<< nameFromAccessionNumber(idmap[iter->entry.id])<<endl;;
    //string found = doStitch(fhpos, iter->second.nucs.substr(0, 100), iter->second.nucs.substr(iter->second.nucs.length()-100), 1.2*iter->second.nucs.length(), 35, false);
    //cout <<"Diff: "<<dnaDiff(found, iter->second.nucs)<<endl;
  }
}
