#pragma once
#include <string>
#include <stdio.h>
#include <stdint.h>
#include <stdexcept>
#include "zstuff.hh"

//! Represents a FastQRead. Can be reversed or not. 
struct FastQRead
{
  FastQRead() : reversed(false), position(0) {}
  std::string d_nucleotides;
  std::string d_quality;
  std::string d_header;
  bool exceedsQuality(unsigned int);
  std::string getSangerQualityString() const;
  void reverse();
  bool reversed;
  uint64_t position; //!< Position in the source file. The 64 bits may encode the file too, it is not a number for the end user to use. Feed it to a FastQReader.
};

//! Reads a single FASTQ file, and can seek in it. Does adapation of quality scores (Sanger by default) and and can also snip off first n or last n bases.
class FASTQReader
{
public:
  FASTQReader(const std::string& str, unsigned int qoffset, unsigned int snipLeft=0, unsigned int snipRight=0);

  void seek(uint64_t pos) 
  {
    d_reader.seek(pos);
  }

  unsigned int getRead(FastQRead* fq); //!< Get a FastQRead, return number of bytes read
private:
  unsigned int d_qoffset;
  unsigned int d_snipLeft, d_snipRight;
  ZLineReader d_reader;
};

//! Reads FASTQs from two (synchronised) files at a time. Does magic with 64 bits offsets to encode which of the two FASTQReader to read from.
class StereoFASTQReader
{
public:
  StereoFASTQReader(const std::string& name1, const std::string& name2, 
		    unsigned int qoffset, unsigned int snipLeft=0, unsigned int snipRight=0) : d_fq1(name1, qoffset, snipLeft, snipRight), d_fq2(name2, qoffset, snipLeft, snipRight) 
  {}

  void seek(uint64_t pos);

  unsigned int getRead(uint64_t pos, FastQRead* fq2);
  unsigned int getReadPair(FastQRead* fq1, FastQRead* fq2);
private:
  FASTQReader d_fq1, d_fq2;
  constexpr static uint64_t s_mask = ~(1ULL<<63);
};
