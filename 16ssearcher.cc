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
#include "githash.h"
#include "dnamisc.hh"
#include <tclap/CmdLine.h>
#include <boost/lexical_cast.hpp>
#include <fstream>
#include <numeric>
#include <set>

extern "C" {
#include "hash.h"

}
#include "misc.hh"
using namespace std;

static map<string, string> g_namescache;

string nameFromAccessionNumber(const std::string& number)
{
  if(g_namescache.count(number))
    return g_namescache[number];
  
  FILE* fp = popen((string("wget -q -O- 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=docsum&tool=antonie&id=")+number+"' | grep Title | cut -f2 -d'>' | cut -f1 -d'<'").c_str(), "r");
  
  if(!fp) 
    throw runtime_error("Unable to open wget pipe: "+string(strerror(errno)));

  string line;
  if(stringfgets(fp, &line)) {
    boost::trim_right(line);
    g_namescache[number]=line;
  }
  pclose(fp);
  return line;
}


Search16S::Search16S(const std::string& fname)
{
  d_linereader = LineReader::make(fname);
}

bool Search16S::get(Entry* entry)
{
  entry->nucs.clear();
  char line[16384];
 
  if(!d_linereader->fgets(line, sizeof(line)))
    return false;

  if(line[0]!='>') 
    throw runtime_error("Unable to parse line '"+string(line)+"' as green genes 16s line, should have >");

  if(line[1]=='S') {
    entry->id = atoi(line+2);
    auto begin = strchr(line, ' '), end = strchr(line, '\t');
    if(begin && end) {
      *begin=0;
      entry->name.assign(begin+1, end);
    }
  }
  else
    entry->id = atoi(line+1);

  for(;;) {
    if(!d_linereader->fgets(line, sizeof(line)))
      return false;
    if(line[0]=='>') {
      d_linereader->unget(line);
      break;
    }

    auto len=strlen(line);
    line[len-1]=0;
    for(unsigned int n=0; n< len-1;++n)
      line[n]=toupper(line[n]);

    char c=line[0];
    if(c<'A' || c>'Z')  // K??
      throw runtime_error("Unable to parse line '"+string(line)+"' as green genes 16s line");

    entry->nucs.append(line);
  }

  return true;
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
  TCLAP::CmdLine cmd("Command description message", ' ', "g" + string(g_gitHash));

  TCLAP::ValueArg<std::string> mode("m","mode","Database mode, gg for Green Genes, rdp for Ribosome Database Project",true,"","mode", cmd);
  TCLAP::ValueArg<int> qualityOffsetArg("q","quality-offset","Quality offset in fastq. 33 for Sanger.",false, 33,"offset", cmd);
  //  TCLAP::ValueArg<int> beginSnipArg("b","begin-snip","Number of nucleotides to snip from begin of reads",false, 0,"nucleotides", cmd);
  //  TCLAP::ValueArg<int> endSnipArg("e","end-snip","Number of nucleotides to snip from end of reads",false, 0,"nucleotides", cmd);
  TCLAP::ValueArg<int> gaps("g","gaps","Number of unmatched bases in 16S fragment allowed",false, 0,"nucleotides", cmd);

  TCLAP::UnlabeledMultiArg<string> multi("filenames", "FASTQ filenames", true, "files",  cmd);
  cmd.parse(argc, argv);
  map<int, string> ggmap;
  map<FASTQReader*, unique_ptr<vector<HashedPos> > > fhpos;

  FASTQReader* fqreader;

  vector<string> files = multi.getValue();
  auto iter = files.begin();
  Search16S s16(*iter++);
  if(mode.getValue()=="gg")
    ggmap=GreenGenesToGenbank(*iter++);
  else if(mode.getValue()!="rdp") {
    cerr<<"Mode needs to be 'gg' or 'rdp'!"<<endl;
    return EXIT_FAILURE;
  }
  for(;iter != files.end(); ++iter) {
    fqreader = new FASTQReader(*iter, qualityOffsetArg.getValue());
    fhpos[fqreader]=indexFASTQ(fqreader, *iter, 35);
  }

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
      //if(candidate.entry.id==128433 || candidate.entry.id==2503772)
      //  continue;
      if(n > 20 && !accumulate(qscores.begin()+n-20, qscores.begin()+n, 0))
	break;
    }

    auto sum=accumulate(qscores.begin(), qscores.end(), 0);
    int holes=0;
    for(unsigned int qpos = 0; qpos < qscores.size(); ++qpos) {
      if(qpos > 20 && !qscores[qpos])
	holes++;
    }
    if(holes <= gaps.getValue()) {
      ofstream cov(boost::lexical_cast<string>(candidate.entry.id)+".cov");
      for(unsigned int qpos = 0; qpos < qscores.size(); ++qpos) {
	cov<<qpos<<'\t'<< qscores[qpos]<<'\n';
      }
      
      cerr<<"\nGood overage ("<<sum/candidate.entry.nucs.size()<<", holes="<<holes<<") on " <<
	candidate.entry.id<<" -> "<<ggmap[candidate.entry.id]<< ": "<< (candidate.entry.name.empty() ? nameFromAccessionNumber(ggmap[candidate.entry.id]) : candidate.entry.name) <<endl;
      candidate.score=sum/candidate.entry.nucs.size();

      candidates.push_back(candidate);
	  
      if(sum >= maxscore) {
	maxscore = sum;
	cerr<<" -> Best current guess: "<<candidate.entry.id<<" ("<<sum/candidate.entry.nucs.size()<<") -> "<<ggmap[candidate.entry.id]<<endl;
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
    cout<<"Score: "<<iter->score<<" ("<<iter->reads.size()*100.0/allReads.size()<<"% of reads), db ID: "<<iter->entry.id<< " -> "<<ggmap[iter->entry.id]<<": "<< (iter->entry.name.empty() ? nameFromAccessionNumber(ggmap[iter->entry.id]) : iter->entry.name)<<endl;
    //string found = doStitch(fhpos, iter->second.nucs.substr(0, 100), iter->second.nucs.substr(iter->second.nucs.length()-100), 1.2*iter->second.nucs.length(), 35, false);
    //cout <<"Diff: "<<dnaDiff(found, iter->second.nucs)<<endl;
  }
}
