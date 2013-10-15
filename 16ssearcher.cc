#include "16ssearcher.hh"
#include <zlib.h>
#include <stdexcept>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <iostream>
#include <boost/foreach.hpp>
#include "fastq.hh"

extern "C" {
#include "hash.h"

}
#include "misc.hh"
using namespace std;

string compactify(const std::string& text)
{
  int pos;
  unsigned char c=0;
  unsigned char t=0;
  string ret;
  string::size_type i;
  for(i = 0; i < text.length(); ++i) {
    if(text[i]=='A') 
      t=0;
    else if(text[i]=='C')
      t=1;
    else if(text[i]=='G')
      t=2;
    else if(text[i]=='T')
      t=3;
    c<<=2;
    c|=t;
    if(i && !(i%4)) {
      ret.append(1,c);
      c=0;
    }
  }
  
  return ret;
}

#define INDEXSIZE 24

Search16S::Search16S(const std::string& fname, int indexlength)
{
  gzFile fp = gzopen(fname.c_str(), "r");
  if(!fp) 
    throw runtime_error("Unable to open file '"+fname+"' for reading");

  char line[16384];
  int mode=0;
  Entry16S e16s;
  unsigned int lineno=0, hashes=0;
  while(gzgets(fp, line, sizeof(line))) {
    if(!((lineno++)%1000)) {
      fprintf(stderr, "\r%d", lineno);
    }
    if((!mode && line[0]!='>') || 
       (mode && line[0]!='C' && line[0]!='A' && line[0]!='G' && line[0]!='T' && line[0]!='K' && line[0]!='N' && line[0]!='M' && line[0]!='B'))  // K??
      throw runtime_error("Unable to parse line '"+string(line)+"' as green genes 16s line");

    if(!mode) {
      e16s.id = atoi(line+1);
      mode=1;
    }
    else {
      mode=0;
      e16s.hits=0;
      chomp(line);
      unsigned int linelen=strlen(line);
      if(linelen < INDEXSIZE) {
	cerr<<"Unable to index short line of length "<<linelen<<endl;
	continue;
      }
      
      for(string::size_type i = 0; i < linelen - INDEXSIZE; i+=80) {
	hashes++;
	uint32_t thehash=hash(line+i, INDEXSIZE, 0);
	d_hashes[thehash].push_back(d_entries.size());
      }
      //      cout<<e16s.id<<" inserted at offset "<<d_entries.size()<<endl;
      d_entries.push_back(e16s);
    }
  }
  cerr<<endl<<"Indexed "<<d_entries.size()<<" 16S entries, resulted in "<<d_hashes.size()<<" different hashes, average fill: "<<hashes/d_hashes.size()<<endl;
}

void Search16S::score(const std::string& nucleotides)
{
  typedef map<uint32_t, uint32_t> hits_t;
  hits_t hits;
  for(string::size_type i = 0; i < nucleotides.size()-INDEXSIZE; ++i) {
    uint32_t hashval = hash(nucleotides.c_str()+i, INDEXSIZE, 0);
    hashes_t::const_iterator iter = d_hashes.find(hashval);
    if(iter == d_hashes.end())
      continue;
    BOOST_FOREACH(unsigned int i, iter->second) {
      hits[i]++;
    }
    BOOST_FOREACH(const hits_t::value_type & val, hits) {
      if(val.second == 2) {
	d_entries[val.first].hits++;
	//	cout<<val.second<<" matches on "<<d_entries[val.first].id<<" ["<<val.first<<"], hits now: "<<d_entries[val.first].hits<<endl;
      }
    }
  }
}

vector<Entry16S> Search16S::topScores()
{
  typedef multimap<uint32_t, const Entry16S*> scores_t;
  scores_t scores;
  BOOST_FOREACH(const Entry16S& e16s, d_entries) {
    scores.insert(make_pair(e16s.hits, &e16s));
  }
  vector<Entry16S> ret;
  for(scores_t::reverse_iterator iter = scores.rbegin(); iter != scores.rend(); ++iter) {
    ret.push_back(*iter->second);
  }
  return ret;
}

void Search16S::printTop(unsigned int limit) 
{
  vector<Entry16S> top = topScores();

  BOOST_FOREACH(const Entry16S& e16s, top) {
    cout<<e16s.id<<": "<<e16s.hits<<endl;
    if(!limit--)
      break;
  }
  cout<<"---"<<endl;

}

int main(int argc, char** argv)
{
  Search16S s16("./gg_13_5.fasta.gz", 150);
  cout<<"Done reading"<<endl;
  FASTQReader fq("tot.fastq");
  FastQRead fqfrag;
  unsigned int num=0;
  while(fq.getRead(&fqfrag)) {
    s16.score(fqfrag.d_nucleotides);
    fqfrag.reverse();
    s16.score(fqfrag.d_nucleotides);
    if(!((num++) % 5000)) {
      s16.printTop(20);
    }
  }


}
