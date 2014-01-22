#include "refgenome.hh"
#include "geneannotated.hh"
#include <iostream>
#include "misc.hh"
#include <set>
#include <algorithm>
#include "dnamisc.hh"
#include <boost/lexical_cast.hpp>
#include <fstream>

extern "C" {
#include "hash.h"
}
using namespace std;

struct HashedPos
{
  uint32_t hash;
  uint64_t position;
  bool operator<(const HashedPos& b) const {
    return hash < b.hash;
  } 
  bool operator<(uint32_t h) const {
    return hash < h;
  }
}__attribute__((packed));

int g_chunklen=35;

unique_ptr<vector<HashedPos> > indexFASTQ(FASTQReader* fqreader, const std::string& fname)
{

  unique_ptr<vector<HashedPos> > hpos(new vector<HashedPos>());
  FILE* fp=fopen((fname+".index").c_str(), "r");
  
  if(fp) {
    auto size=filesize((fname+".index").c_str());
    if(size % sizeof(HashedPos)) {
      fclose(fp);
      throw runtime_error("Index has wrong size. Sizeof(HashedPos): "+boost::lexical_cast<string>(sizeof(HashedPos)));
    }
    hpos->resize(size/sizeof(HashedPos));
    if(fread(&(*hpos)[0], 1, size, fp) != size)
      throw runtime_error("Index corrupt");
    fclose(fp);
    return hpos;
  }
  cerr<<"Indexing "<<fname<<endl;
  FastQRead fqr;
  while(fqreader->getRead(&fqr)) {
    uint32_t h = hash(fqr.d_nucleotides.c_str(), g_chunklen, 0);
    hpos->push_back({h, fqr.position});
    fqr.reverse();
    h = hash(fqr.d_nucleotides.c_str(), g_chunklen, 0);
    hpos->push_back({h, fqr.position});
  }
  std::sort(hpos->begin(), hpos->end());

  fp=fopen((fname+".index").c_str(), "w");
  for(const auto& hpo : *hpos) {
    fwrite(&hpo.hash, 1, sizeof(hpo.hash), fp);
    fwrite(&hpo.position, 1, sizeof(hpo.position), fp);
  }
  fclose(fp);

  return hpos;
}
int g_maxdepth;

dnapos_t g_record=0;

set<pair<uint32_t, dnapos_t> > g_beenthere;

string g_bestcontig;

set<string> g_candidates;

int dnaDiff(const std::string& a, const std::string& b)
{
  if(a==b)
    return 0;
  //  cout<<"A: "<<a<<endl;
  //  cout<<"B: "<<b<<endl;
  //  cout<<"   ";
  int diff=0;
  for(string::size_type o = 0 ; o < a.length() && o < b.length(); ++o) {
    if(a[o] != b[o]) {
      //      cout<<'!';
      diff++;
    }
    else
      ;
      //      cout<<' ';
  }
  
  //  cout<<endl<<"diff: "<<diff<<endl;
  return diff;
}

vector<FastQRead> getConsensusMatches(const std::string& consensus, map<FASTQReader*, unique_ptr<vector<HashedPos> > >& fhpos)
{
  vector<FastQRead> ret;
  if(consensus.find('N') != string::npos)
    return ret;

  uint32_t h = hash(consensus.c_str(), consensus.length(), 0);
  //  cout<<"Looking for "<<consensus<<endl;
  HashedPos fnd({h, 0});
  
  vector<FastQRead> options;
  
  for(auto& hpos : fhpos) {
    auto range = equal_range(hpos.second->begin(), hpos.second->end(), fnd);
    for(;range.first != range.second; ++range.first) {
      FastQRead fqr;
      //cout<<"\tFound potential hit at offset "<<range.first->position<<"!"<<endl;
      hpos.first->seek(range.first->position);
      hpos.first->getRead(&fqr);
      if(fqr.d_nucleotides.substr(0,g_chunklen) != consensus) {
	fqr.reverse();
	
	if(fqr.d_nucleotides.substr(0,g_chunklen) != consensus) {
	  cout<<"\thash lied\n";
	  continue;
	}
      }
      ret.push_back(fqr);
    }
  }
  return ret;
}

struct Base
{
  Base() : aCount(0), cCount(0), gCount(0), tCount(0){};
  int aCount, cCount, gCount, tCount;
  char getBest()
  {
    if(!aCount && !cCount && !gCount && !tCount)
      return 'N';
    vector<pair<int, int*> > scores{
      {aCount, &aCount},
	{cCount, &cCount},
	  {gCount, &gCount},
	    {tCount, &tCount}};
    sort(scores.begin(), scores.end());
    auto& best = scores[3].second;
    if(best == &aCount)
      return 'A';
    else if(best == &cCount)
      return 'C';
    else if(best == &gCount)
      return 'G';
    else 
      return 'T';
  }
  void feed(char c, int amount=1)
  {
    if(c=='A')
      aCount+=amount;
    else if(c=='C')
      cCount+=amount;
    else if(c=='G')
      gCount+=amount;
    else if(c=='T')
      tCount+=amount;
  }
};


// stitcher fasta startpos fastq fastq
int main(int argc, char**argv)
{
  ReferenceGenome rg(argv[1]);

  string startseed;
  if(isalpha(argv[2][0]))
    startseed = argv[2];
  else {
    dnapos_t startpos = atoi(argv[2]);
    startseed = rg.snippet(startpos, startpos+100);
  }

  dnapos_t endpos = atoi(argv[3]);
  map<FASTQReader*, unique_ptr<vector<HashedPos> > > fhpos;

  FASTQReader* fqreader;

  for(int f = 4; f < argc; ++f) {
    fqreader = new FASTQReader(argv[f], 33, 0);
    fhpos[fqreader]=indexFASTQ(fqreader, argv[f]);
  }
  string genome;

  string endseed = rg.snippet(endpos, endpos+100);

  cout << "Startseed: "<<startseed<< " (" <<startseed.size()<<")\n";
  cout << "Endseed: "<<endseed<<endl;
  //  cout << "Reference: \n"<<rg.snippet(startpos, endpos+100) << endl;

  cout << startseed<<endl;
  int offset=0;
  string totconsensus;
  // cons:
  // start: ABCDEFGHIJKLMNOPQRSTUVWXYZ
  // new:          HIJKLMNOPQRSTYVWXYZ123456
  // new:                 OPQRSTYVWXYZ123456789
  // new:                    RSTYVWXYZ123456789012
  // new                          WXYZ123456789012ABCDEF
  // end:   ABCDEFGHIJKLMNOPQRSTUVWXYZ123456789012  // <- 50% longer
  // cons:  ABCDEFGHIJKLMNOPQRSTUVWXYZ // move 100% of original length of consensus so connsensus
  //
  // start: NOPQRSTUVWXYZ123456789012  // move ahead 50% of original lenght
  for(;;) {
    vector<pair<string,string> > story;
    story.push_back(make_pair(startseed, string(startseed.size(), (char)40)));

    for(unsigned int n=0; n < startseed.size() - g_chunklen;++n) {
      string part=startseed.substr(n, g_chunklen);
      auto matches = getConsensusMatches(part, fhpos);
      for(auto& match : matches) {
	int diff = dnaDiff(startseed.substr(n), match.d_nucleotides);
	if(diff < 5) {
	  cout << string(offset,'-')<<string(n, ' ') << match.d_nucleotides<<endl;
	  story.push_back({string(n, ' ')+match.d_nucleotides,
		string(n, ' ')+match.d_quality});
	}
      }
    }
    vector<Base> consensus;
    consensus.resize(startseed.size()*1.5);
    for(unsigned int n = 0 ; n < consensus.size(); ++n) {
      for(const auto& candidate : story) {
	if(n < candidate.first.size())
	  consensus[n].feed(candidate.first[n], candidate.second[n]);
      }
    }
    string newconsensus;
    cout<<totconsensus;
    for(unsigned int n = 0 ; n < consensus.size();++n) {
      cout<<consensus[n].getBest();
      newconsensus.append(1, consensus[n].getBest());
    }
    cout<<endl;
    startseed=newconsensus.substr(startseed.length()/2, startseed.length());
    totconsensus+=newconsensus.substr(0, startseed.length()/2);
    offset+=startseed.length()/2;
    cout<<"--"<<endl;
  }

  
}

