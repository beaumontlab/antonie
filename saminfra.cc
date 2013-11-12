#include "saminfra.hh"
#include "fastq.hh"
#include <stdexcept>
#include <string.h>
#include <string>
#include <boost/lexical_cast.hpp>
using std::string;
using boost::lexical_cast;

SAMWriter::~SAMWriter()
{
  if(d_fp)
    fclose(d_fp);
}

SAMWriter::SAMWriter(const std::string& fname, const std::string& genomeName, dnapos_t len) : d_fname(fname), d_genomeName(genomeName)
{
  if(fname.empty()) {
    d_fp=0;
    return;
  }
  d_fp = fopen(fname.c_str(), "w");
  if(!d_fp) 
    throw std::runtime_error("Unable to open '"+fname+"' for writing SAM file"+strerror(errno));

  fprintf(d_fp, "@HD\tVN:1.0\tSO:unsorted\n");
  fprintf(d_fp, "@SQ\tSN:%s\tLN:%u\n", d_genomeName.c_str(), len);
  fprintf(d_fp, "@PG\tID:antonie\tPN:antonie\tVN:0.0.0\n");  
}

void SAMWriter::write(dnapos_t pos, const FastQRead& fqfrag, int indel)
{
  if(!d_fp) 
    return;

  string header;
  string::size_type spacepos = fqfrag.d_header.find(' ');
  if(spacepos != string::npos)
    header = fqfrag.d_header.substr(0, spacepos);
  else
    header = fqfrag.d_header;

  string quality = fqfrag.d_quality;
  for(auto& c : quality) {
    c+=33; // we always output Sanger
  }

  string cigar;
  if(!indel) {
    cigar = lexical_cast<string>(fqfrag.d_nucleotides.length());
    cigar.append(1,'M');
  }
  else if(indel < 0) {
    cigar = lexical_cast<string>(-indel);
    cigar.append(1,'M');
    cigar += "1D";
    cigar += lexical_cast<string>(fqfrag.d_nucleotides.length()+indel);
    cigar.append(1,'M');
  }
  else if(indel > 0) {
    cigar = lexical_cast<string>(indel);
    cigar.append(1,'M');
    cigar += "1I";
    cigar += lexical_cast<string>(fqfrag.d_nucleotides.length()-1-indel);
    cigar.append(1,'M');

  }
	
  fprintf(d_fp, "%s\t%u\t%s\t%u\t42\t%s\t*\t0\t0\t%s\t%s\n",
	  header.c_str(), 
	  fqfrag.reversed ? 0x10: 0,
	  d_genomeName.c_str(), pos, cigar.c_str(),
	  fqfrag.d_nucleotides.c_str(), quality.c_str());  
}


