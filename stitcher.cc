#include "refgenome.hh"
#include "geneannotated.hh"
#include <iostream>
#include "misc.hh"
#include <set>
#include <algorithm>
#include "dnamisc.hh"
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
};

int g_chunklen=30;

unique_ptr<vector<HashedPos> > indexFASTQ(FASTQReader* fqreader, const std::string& fname)
{

  unique_ptr<vector<HashedPos> > hpos(new vector<HashedPos>());
  FILE* fp=fopen((fname+".index").c_str(), "r");
  if(fp) {
    HashedPos hpo;
    while(fread(&hpo.hash, 1, sizeof(hpo.hash), fp) == sizeof(hpo.hash) &&  
	  fread(&hpo.position, 1, sizeof(hpo.position), fp) == sizeof(hpo.position))
      hpos->push_back(hpo);
    if(!feof(fp)) 
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

struct Hypothesis
{
  Hypothesis(const std::string& str, Hypothesis* parent=0) : d_nucs(str), d_trypos(20), d_parent(parent)
  {}
  string d_nucs; 
  unsigned int d_trypos;

  vector<pair<int, Hypothesis*> > d_children;
  Hypothesis* d_parent;
  
  bool isDead()
  {
    return d_trypos+g_chunklen >= d_nucs.length();
  }
  void getLiveLeaves(vector<Hypothesis*>& ret);
  string tracePrint();
  pair<int, Hypothesis*> getLongestLiveLeave(dnapos_t len=0);
  int depth(int d=0) const;
  int print(int d=0) const;
  int longest(int d=0) const;
};


string Hypothesis::tracePrint()
{
  int ourpos=-1;
  auto parent = d_parent;
  string ret = d_nucs;
  vector<pair<int, Hypothesis*> > chain;
  Hypothesis* us = this; 
  while(parent) {
    for(const auto& child : parent->d_children) {
      if(child.second == us) {
	ourpos = child.first;
        chain.push_back({ourpos, parent});
	break;
      }
    }
    ret = parent->d_nucs.substr(0, ourpos) + " " + ret;
    us = parent;
    parent = parent->d_parent; 
  }
  /*
  reverse(chain.begin(), chain.end());
  int offset=0;
  for(auto& e : chain) {
    cout<<string(offset, ' ')<<e.second->d_nucs<<endl;
    offset+=e.first;
  }
  */
  return ret;

} 

int Hypothesis::depth(int d) const
{
  set<int> depths;
  d++;
  depths.insert(d);
  for(const auto& child : d_children) 
    depths.insert(child.second->depth(d));
  return *depths.rbegin();
}

pair<int, Hypothesis*> Hypothesis::getLongestLiveLeave(dnapos_t len) 
{
  set<pair<int, Hypothesis*> > depths;
  depths.insert({len, this});
  for(const auto& child : d_children) 
    depths.insert(child.second->getLongestLiveLeave(len + child.first));
  for(auto iter = depths.rbegin(); iter != depths.rend(); ++iter) 
    if(iter->second && !iter->second->isDead())
      return *iter;
  return {0,nullptr};
}


int Hypothesis::print(int offset) const
{
  cout<<string(offset, ' ')<<d_nucs<<endl;
  for(const auto& child : d_children) 
    child.second->print(offset+child.first);
  return 0;
}

int Hypothesis::longest(int offset) const
{
  int maxret=offset;
  for(const auto& child : d_children) {
    maxret = max(child.second->longest(offset+child.first), maxret);
  }
  return maxret;
}

void Hypothesis::getLiveLeaves(vector<Hypothesis*>& ret)
{
  if(d_trypos != d_nucs.length())
    ret.push_back(this);
  for(auto& child : d_children)
    child.second->getLiveLeaves(ret);
}

vector<string> getConsensusMatches(const std::string& consensus, map<FASTQReader*, unique_ptr<vector<HashedPos> > >& fhpos)
{
  vector<string> ret;
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
      ret.push_back(fqr.d_nucleotides);
    }
  }
  return ret;
}

// stitcher fasta startpos fastq fastq
int main(int argc, char**argv)
{
  ReferenceGenome rg(argv[1]);
  dnapos_t startpos = atoi(argv[2]);
  dnapos_t endpos = atoi(argv[3]);
  map<FASTQReader*, unique_ptr<vector<HashedPos> > > fhpos;

  FASTQReader* fqreader;

  for(int f = 4; f < argc; ++f) {
    fqreader = new FASTQReader(argv[f], 33, 10);
    fhpos[fqreader]=indexFASTQ(fqreader, argv[f]);
  }
  string genome;
  string startseed = rg.snippet(startpos, startpos+100);
  string endseed = rg.snippet(endpos, endpos+100);

  cout << "Startseed: "<<startseed<<endl;
  cout << "Endseed: "<<endseed<<endl;
  cout << "Reference: "<<rg.snippet(startpos, endpos+100) << endl;

  Hypothesis hypo(startseed);
  for(unsigned int n=0;;++n) {
    auto longest=hypo.getLongestLiveLeave();
    
    Hypothesis* leaf = longest.second;
    if(!leaf)
      break;
    
    if(!(n%100))
      cout<<"Interim: "<<leaf->tracePrint()<<endl;
    
    string consensus(leaf->d_nucs.substr(leaf->d_trypos, g_chunklen));

    auto matches = getConsensusMatches(consensus, fhpos);
    for(auto& match : matches) {
      int diff = dnaDiff(leaf->d_nucs.substr(leaf->d_trypos), match);
      if(diff < 5) {
	auto newChild = new Hypothesis(match, leaf);
	leaf->d_children.push_back({leaf->d_trypos, newChild});	  
	if(match.find(endseed) != string::npos) {
	  cout<<"Got hit: "<<endl;
	  cout<<newChild->tracePrint()<<endl;
	  leaf->d_trypos=leaf->d_nucs.length();
	  newChild->d_trypos=newChild->d_nucs.length();
	}
      }
    }
    leaf->d_trypos++;
    
    cout<<"Longest: "<<hypo.longest()<<endl;
    //    if(!(n%100))
    //  hypo.print();
  }
  
  cout<<"Candidates: "<<endl;
  for(auto& candidate : g_candidates) {
    cout<<candidate<<endl;
  }
  rg.index(g_chunklen);
  
  ofstream graph("synten");
  for(auto iter = g_bestcontig.begin(); iter + g_chunklen < g_bestcontig.end(); ++iter) {
    string stretch(iter, iter+g_chunklen);
    for(int n=0; n < 2; ++n) {
      auto positions = rg.getReadPositions(stretch);
      for(auto pos : positions) {
	graph << (startpos + (iter - g_bestcontig.begin())) << '\t';
	if(!n)
	  graph<< "NaN\t"<<pos << '\n';
	else
	  graph << pos << "\tNaN\n";
      }
      
      reverseNucleotides(&stretch);
    }
  }
}



#if 0
  while(genome.length() < rg.size()) {
    startProcessing(snippet, 0, endseed, fhpos);
    g_bestcontig = g_bestcontig.substr(0, g_bestcontig.size()-1000);
    genome += g_bestcontig;
    newgenome << g_bestcontig;

    snippet = genome.substr(genome.length()-100);
    g_record=0;
    g_beenthere.clear();
    g_bestcontig.clear();
    cerr<<"genome size: "<<genome.size()<<endl;
    cerr<<"New snippet: "<<snippet<<endl;
  }
#endif

#if 0
struct Hyposition
{
  Hyposition() : aCount(0), cCount(0), gCount(0), tCount(0)
  {}
  int aCount, cCount, gCount, tCount;
};

struct Hypothesis
{
  void mapString(const std::string& str, dnapos_t position);
  string getConsensus(dnapos_t position, dnapos_t length);
  string getConsensus() 
  {
    return getConsensus(0, d_hypo.size());
  }
  void getEnsemble(dnapos_t position, dnapos_t length, vector<string>* ensemble, string sofar="");
  int diff(dnapos_t position, const std::string& sug);
  vector<Hyposition> d_hypo;
  size_t size()
  {
    return d_hypo.size();
  }
};

void Hypothesis::mapString(const std::string& str, dnapos_t position)
{
  if(d_hypo.size() < position+str.length())
    d_hypo.resize(position+str.length());
  for(string::size_type i = 0; i < str.length(); ++position, ++i) {
    Hyposition& hyp = d_hypo[position];
    acgtDo(str[i],
	       [&]() { hyp.aCount++; }, 
	       [&]() { hyp.cCount++; }, 
	       [&]() { hyp.gCount++; }, 
	       [&]() { hyp.tCount++; });

  }
}

string Hypothesis::getConsensus(dnapos_t position, dnapos_t length)
{
  string ret;
  dnapos_t endpos = position+length;
  for(; position  < endpos; ++position) {
    Hyposition& hyp = d_hypo[position];
    vector<pair<int, char>> counts{{hyp.aCount, 'A'}, 
	{hyp.cCount, 'C'}, 
	  {hyp.gCount, 'G'}, 
	    {hyp.tCount, 'T'} };
  
    sort(counts.begin(), counts.end());
    ret.append(1, counts.rbegin()->second);
    for(auto c: counts)
      cout<< c.second<<": "<<c.first<< " ";
    cout<<endl;
  }  
  return ret;
}

void Hypothesis::getEnsemble(dnapos_t position, dnapos_t length, vector<string>* ensemble, string sofar)
{
  if(ensemble->size() > 5000000)
    return;
  //  cout<<"Position is position "<<position<<" length is "<<length<<endl;
  string ret;
  dnapos_t endpos = position+length;
  for(; position  < endpos; ++position) {
    Hyposition& hyp = d_hypo[position];
    if(hyp.aCount && !hyp.cCount && !hyp.gCount && !hyp.tCount)
      sofar+="A";
    else
    if(!hyp.aCount && hyp.cCount && !hyp.gCount && !hyp.tCount)
      sofar+="C";
    else
    if(!hyp.aCount && !hyp.cCount && hyp.gCount && !hyp.tCount)
      sofar+="G";
    else
    if(!hyp.aCount && !hyp.cCount && !hyp.gCount && hyp.tCount)
      sofar+="T";
    else {
      if(hyp.aCount) {
	getEnsemble(position+1, endpos - position - 1, ensemble, sofar+"A");
      }
      if(hyp.cCount) {
	getEnsemble(position+1, endpos - position - 1, ensemble, sofar+"C");
      }
      if(hyp.gCount) {
	getEnsemble(position+1, endpos - position - 1, ensemble, sofar+"G");
      }
      if(hyp.tCount) {
	getEnsemble(position+1, endpos - position - 1, ensemble, sofar+"T");
      }
      return;
    }    
  }
  ensemble->push_back(sofar);

}

int Hypothesis::diff(dnapos_t position, const std::string& sug)
{
  dnapos_t endpos = position + sug.length();
  if(endpos > d_hypo.size())
    endpos = d_hypo.size();

  int diff=0;
  for(string::size_type sugpos=0; sugpos < sug.length() && position < endpos ; ++position, ++sugpos) {
    Hyposition& hyp = d_hypo[position];
    vector<pair<int, char>> counts{{hyp.aCount, 'A'}, 
	{hyp.cCount, 'C'}, 
	  {hyp.gCount, 'G'}, 
	    {hyp.tCount, 'T'} };
  
    sort(counts.begin(), counts.end());
    if(sug[sugpos] != counts[3].second && !(counts[2].first && sug[sugpos] == counts[2].second))
      diff++;
  }
  return diff;
}
#endif

#if 0
bool startProcessing(Hypothesis& hypo, dnapos_t pos, const std::string& endseed, const map<FASTQReader*, unique_ptr<vector<HashedPos> > > & fhpos, int depth=0)
{
  bool reallyDone=false;
  for( ; pos < hypo.size()-g_chunklen; ++pos) {
    vector<string> ensemble;
    
    hypo.getEnsemble(pos, g_chunklen, &ensemble);
    cout<<"Got "<<ensemble.size()<<" ensemble options"<<endl;
    for(auto& consensus : ensemble) {
      uint32_t h = hash(consensus.c_str(), consensus.length(), 0);
      cout<<"At position "<<pos<<", looking for "<<consensus<<endl;
      HashedPos fnd({h, 0});
      
      vector<FastQRead> options;
      
      for(auto& hpos : fhpos) {
	auto range = equal_range(hpos.second->begin(), hpos.second->end(), fnd);
	for(;range.first != range.second; ++range.first) {
	  FastQRead fqr;
	  cout<<"\tFound potential hit at offset "<<range.first->position<<"!"<<endl;
	  hpos.first->seek(range.first->position);
	  hpos.first->getRead(&fqr);
	  if(fqr.d_nucleotides.substr(0,g_chunklen) != consensus) {
	    fqr.reverse();
	    
	    if(fqr.d_nucleotides.substr(0,g_chunklen) != consensus) {
	      cout<<"\thash lied\n";
	      continue;
	    }
	  }
	  int diff = hypo.diff(pos, fqr.d_nucleotides);
	  if(diff < 15) {
	    cout<<"Mapping @ "<<pos<<", diff = "<<diff<<endl;
	    hypo.mapString(fqr.d_nucleotides, pos);
	  }
	  else
	    cout<<"Not mapping, diff = "<<diff<<endl;
	  
	}
      }
    }
  }
  return reallyDone;
}
#endif
