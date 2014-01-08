#include "refgenome.hh"
#include "geneannotated.hh"
#include <iostream>
#include "misc.hh"
#include <set>
#include <algorithm>
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
    uint32_t h = hash(fqr.d_nucleotides.c_str(), 50, 0);
    hpos->push_back({h, fqr.position});
    fqr.reverse();
    h = hash(fqr.d_nucleotides.c_str(), 50, 0);
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
  cout<<"A: "<<a<<endl;
  cout<<"B: "<<b<<endl;
  cout<<"   ";
  int diff=0;
  for(string::size_type o = 0 ; o < a.length() && o < b.length(); ++o) {
    if(a[o] != b[o]) {
      cout<<'!';
      diff++;
    }
    else
      cout<<' ';
  }
  
  cout<<endl<<"diff: "<<diff<<endl;
  return diff;
}

bool startProcessing(string snippet, dnapos_t pos, const std::string& endseed, const map<FASTQReader*, unique_ptr<vector<HashedPos> > > & fhpos, int depth=0)
{
  if(!depth)
    g_bestcontig.clear();
  bool reallyDone=false;
  if(!(g_maxdepth < depth)) {
  //    cerr<<"New depth: "<<depth<<endl;
    g_maxdepth = depth;
  } 

  string::size_type endseedpos=snippet.find(endseed);
  if(endseedpos != string::npos) {
    g_candidates.insert({snippet.substr(0, endseedpos + endseed.length())});
    cerr<<"Got one, now have "<<g_candidates.size()<<" candidates"<<endl;
    cout<<"Got one, now have "<<g_candidates.size()<<" candidates"<<endl;
    cout<<"Candidate: "<<snippet.substr(0, endseedpos + endseed.length())<<endl;
    return false;
  }

  if(snippet.length() > 3000 ) 
    return true;
  
  for( ; pos < snippet.length()-50; ++pos) {
    uint32_t h = hash(snippet.c_str(), snippet.length(), 0);
    if(g_beenthere.count({h, pos})) {
      cout<<"Been here"<<endl;
      goto done;
    }
    g_beenthere.insert({h,pos});
    
    h = hash(snippet.c_str() + pos, 50, 0);
    cout<<"Looking for hash '"<<h<<"':\n";
    HashedPos fnd({h, 0});

    vector<FastQRead> options;
    
    for(auto& hpos : fhpos) {
      auto range = equal_range(hpos.second->begin(), hpos.second->end(), fnd);
      for(;range.first != range.second; ++range.first) {
	FastQRead fqr;
	cout<<"\tFound potential hit at offset "<<range.first->position<<"!"<<endl;
	hpos.first->seek(range.first->position);
	hpos.first->getRead(&fqr);
	if(fqr.d_nucleotides.substr(0,50) != snippet.substr(pos, 50)) {
	  fqr.reverse();
      
	  if(fqr.d_nucleotides.substr(0,50) != snippet.substr(pos, 50)) {
	    cout<<"\thash lied\n";
	    continue;
	  }
	}
       
	string belief(snippet.substr(pos+50));
	string measure(fqr.d_nucleotides.substr(50, snippet.length()-pos-50));
	auto diff = dnaDiff(belief, measure);
	
	if(diff < 4) {
	  options.push_back(fqr);  
	}	
      }
    }

    set<string> tot;
    cout<<"Have "<<options.size()<<" options at this point"<<endl;
    for(auto& option: options) {
      if(tot.count(option.d_nucleotides)) {
	cout<<" duplicate"<<endl;
	continue;
      }
      tot.insert(option.d_nucleotides);
      auto additlen = pos + option.d_nucleotides.length() - snippet.length();
      cout<<string(pos, ' ')<<option.d_nucleotides<< " ("<<additlen<<", "<<option.reversed<<")"<<endl;
      //      string newbit = option.d_nucleotides.substr(option.d_nucleotides.length() - additlen);
      //cout<<snippet+newbit<<endl;
     
      string newsnippet=snippet.substr(0, pos) + option.d_nucleotides;
      cout<<newsnippet<<endl;
      if(startProcessing(newsnippet, pos+1, endseed, fhpos, depth+1)) {
	reallyDone=true;
	goto done;
      }
    }
  }
 done:;
  if(snippet.length() > g_record) {
    cout<<"Final: "<<snippet.length()<<" "<<snippet<<endl;
    g_record = snippet.length();
    g_bestcontig = snippet;
    cerr<<g_record<<endl;
  }
  return reallyDone;
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

  startProcessing(startseed, 0, endseed, fhpos);

  cout<<"Candidates: "<<endl;
  for(auto& candidate : g_candidates) {
    cout<<candidate<<endl;
  }
  rg.index(50);
  
  ofstream graph("synten");
  for(auto iter = g_bestcontig.begin(); iter + 50 < g_bestcontig.end(); ++iter) {
    string stretch(iter, iter+50);
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
