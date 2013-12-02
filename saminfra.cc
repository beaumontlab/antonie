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

void SAMWriter::write(dnapos_t pos, const FastQRead& fqfrag, int indel, int flags, const std::string& rnext, dnapos_t pnext, int32_t tlen)
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
	
  fprintf(d_fp, "%s\t%u\t%s\t%u\t42\t%s\t"
	  "%s\t%u\t%d\t"
	  "%s\t%s\n",
	  header.c_str(), 
	  flags + (fqfrag.reversed ? 0x10: 0),
	  d_genomeName.c_str(), pos, cigar.c_str(),
	  rnext.c_str(), pnext, tlen,
	  fqfrag.d_nucleotides.c_str(), quality.c_str());  
}


/* calculate bin given an alignment covering [beg,end) (zero-based, half-close-half-open) */
static int reg2bin(int beg, int end)
{
  --end;
  if (beg>>14 == end>>14) return ((1<<15)-1)/7 + (beg>>14);
  if (beg>>17 == end>>17) return ((1<<12)-1)/7 + (beg>>17);
  if (beg>>20 == end>>20) return ((1<<9)-1)/7 + (beg>>20);
  if (beg>>23 == end>>23) return ((1<<6)-1)/7 + (beg>>23);
  if (beg>>26 == end>>26) return ((1<<3)-1)/7 + (beg>>26);
  return 0;
}


struct BAMBuilder
{
  BAMBuilder(std::string* str) : d_str(str)
  {}

  void write(const char* p, int num)
  {
    d_str->append(p, num);
  }

  void write32(uint32_t val) 
  {
    d_str->append((const char*)&val, 4);
  }

  void writeBAMString(const std::string& str)
  {
    write32(str.length());
    d_str->append(str);
  }

  string* d_str;
};

BAMWriter::BAMWriter(const std::string& fname, const std::string& header, const std::string& refname, dnapos_t reflen) 
{
  BAMBuilder bb(&d_block);

  char magic[]="BAM\1";

  bb.write(magic, 4);
  bb.writeBAMString(header);
  bb.write32(1);
  bb.writeBAMString(refname);  
  bb.write32(reflen);
}

