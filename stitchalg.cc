#include "fastqindex.hh"
#include <iostream>
#include <algorithm>
#include "stitchalg.hh"

using namespace std;

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

unsigned int g_chunklen=35;

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



string doStitch(const map<FASTQReader*, unique_ptr<vector<HashedPos> > >& fhpos, const std::string& startseed_,
		const std::string& endseed, unsigned int maxlen, int chunklen)
{
  string startseed(startseed_);
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

    for(unsigned int n=0; n < startseed.size() - chunklen;++n) {
      string part=startseed.substr(n, chunklen);
      auto matches = getConsensusMatches(part, fhpos, chunklen);
      for(auto& match : matches) {
	int diff = dnaDiff(startseed.substr(n), match.d_nucleotides);
	if(diff < 5) {
	  //	  cout << string(offset,'-')<<string(n, ' ') << match.d_nucleotides<<endl;
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
    string::size_type endpos = totconsensus.find(endseed);
    if(endpos != string::npos) {
      totconsensus.resize(endpos+endseed.size());
      cout<<"Done: \n"<<totconsensus<<endl;
      break;
    }
    if(totconsensus.size() > maxlen)
      break;
  }
  return totconsensus;
}
