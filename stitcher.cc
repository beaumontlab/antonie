#include "refgenome.hh"
#include "geneannotated.hh"
#include <iostream>
#include "misc.hh"
#include <set>
#include <algorithm>
#include "dnamisc.hh"
#include <boost/lexical_cast.hpp>
#include "fastqindex.hh"
#include <fstream>
#include "stitchalg.hh"
extern "C" {
#include "hash.h"
}

using namespace std;

int g_maxdepth;

dnapos_t g_record=0;

set<pair<uint32_t, dnapos_t> > g_beenthere;

string g_bestcontig;

set<string> g_candidates;


// stitcher fasta startpos fastq fastq
int main(int argc, char**argv)
{
  ReferenceGenome rg(argv[1]);

  int chunklen=35;
  string startseed;
  if(isalpha(argv[2][0]))
    startseed = argv[2];
  else {
    dnapos_t startpos = atoi(argv[2]);
    startseed = rg.snippet(startpos, startpos+100);
  }

  map<FASTQReader*, unique_ptr<vector<HashedPos> > > fhpos;

  FASTQReader* fqreader;

  for(int f = 4; f < argc; ++f) {
    fqreader = new FASTQReader(argv[f], 33, 0);
    fhpos[fqreader]=indexFASTQ(fqreader, argv[f], chunklen);
  }

  string endseed = argv[3];

  doStitch(fhpos, startseed, endseed, 10000, chunklen, true);
}

