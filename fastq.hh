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
  std::string getNameFromHeader() const;
  bool exceedsQuality(unsigned int);
  std::string getSangerQualityString() const;
  void reverse();
  bool reversed;
  uint64_t position; //!< Position in the source file. The 64 bits may encode the file too, it is not a number for the end user to use. Feed it to a FastQReader.

  bool operator<(const FastQRead& rhs) const
  {
    return std::tie(d_nucleotides, d_quality, reversed, position) < 
      std::tie(rhs.d_nucleotides, rhs.d_quality, rhs.reversed, rhs.position);
  }

};

//! Reads a single FASTQ file, and can seek in it. Does adapation of quality scores (Sanger by default) and and can also snip off first n or last n bases.
class FASTQReader
{
public:
  FASTQReader(const std::string& str, unsigned int qoffset);
  void setTrim(unsigned int trimLeft, unsigned int trimRight)
  {
    d_snipLeft = trimLeft;
    d_snipRight = trimRight;
  }
  void seek(uint64_t pos) 
  {
    d_reader->seek(pos);
  }
  uint64_t estimateReads();
  unsigned int getRead(FastQRead* fq); //!< Get a FastQRead, return number of bytes read
private:
  unsigned int d_qoffset;
  unsigned int d_snipLeft, d_snipRight;
  std::unique_ptr<LineReader> d_reader;
};

//! Reads FASTQs from two (synchronised) files at a time. Does magic with 64 bits offsets to encode which of the two FASTQReader to read from.
class StereoFASTQReader
{
public:
  StereoFASTQReader(const std::string& name1, const std::string& name2, 
		    unsigned int qoffset) : d_fq1(name1, qoffset), d_fq2(name2, qoffset) 
  {}

  void setTrim(unsigned int trimLeft, unsigned int trimRight);
  void seek(uint64_t pos);
  uint64_t estimateReads();
  unsigned int getRead(uint64_t pos, FastQRead* fq2);
  unsigned int getReadPair(FastQRead* fq1, FastQRead* fq2);
private:
  FASTQReader d_fq1, d_fq2;
  static uint64_t s_mask;
};
