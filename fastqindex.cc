#include "fastqindex.hh"
#include "misc.hh"
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <algorithm>
using namespace std;

extern "C" {
#include "hash.h"
}

unique_ptr<vector<HashedPos> > indexFASTQ(FASTQReader* fqreader, const std::string& fname, int chunklen)
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
    uint32_t h = hash(fqr.d_nucleotides.c_str(), chunklen, 0);
    hpos->push_back({h, fqr.position});
    fqr.reverse();
    h = hash(fqr.d_nucleotides.c_str(), chunklen, 0);
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

vector<FastQRead> getConsensusMatches(const std::string& consensus, const map<FASTQReader*, unique_ptr<vector<HashedPos> > >& fhpos, int chunklen)
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
      //      cout<<"\tFound potential hit at offset "<<range.first->position<<"!"<<endl;
      hpos.first->seek(range.first->position);
      hpos.first->getRead(&fqr);
      if(fqr.d_nucleotides.substr(0,chunklen) != consensus) {
	fqr.reverse();
	
	if(fqr.d_nucleotides.substr(0,chunklen) != consensus) {
	  //	  cout<<"\thash lied\n";
	  continue;
	}
      }
      ret.push_back(fqr);
    }
  }
  return ret;
}
