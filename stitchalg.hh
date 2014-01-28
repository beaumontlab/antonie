#pragma once
#include <string>
#include "fastqindex.hh"

std::string doStitch(const std::map<FASTQReader*, std::unique_ptr<std::vector<HashedPos> > >& fhpos, 
		     const std::string& startseed_,
		     const std::string& endseed, unsigned int maxlen, int chunklen, bool verbose);
int dnaDiff(const std::string& a, const std::string& b);
