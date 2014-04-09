#pragma once
#include <string>
#include <zlib.h>
#include <stdio.h>
#include <map>
#include <stdexcept>
#include <boost/utility.hpp>
#include <memory>
#include <boost/crc.hpp>
#include <stdint.h>

//! Virtual base for seekable line readers
class LineReader
{
public:
  virtual ~LineReader() {}
  // virtual bool getLine(std::string* str) = 0;
  virtual char* fgets(char* line, int num) = 0;
  virtual void seek(uint64_t pos) = 0;
  virtual uint64_t getUncPos()=0;
  virtual void unget(char *line) = 0;
  virtual uint64_t uncompressedSize() = 0;
  static std::unique_ptr<LineReader> make(const std::string& fname);
};

//! A plain text seekable line reader
class PlainLineReader : public LineReader, boost::noncopyable
{
public:
  PlainLineReader(const std::string& fname);
  ~PlainLineReader();
  //  bool getLine(std::string* str) = 0;
  char* fgets(char* line, int num);
  void seek(uint64_t pos);
  uint64_t getUncPos();
  void unget(char *line);
  uint64_t uncompressedSize();
private:
  FILE* d_fp;
  std::string d_stash;
};


//! A gzipped compressed seekable line reader
class ZLineReader : public LineReader, boost::noncopyable
{
public:
  ZLineReader(const std::string& fname);
  ~ZLineReader();
  //  bool getLine(std::string* str);
  char* fgets(char* line, int num);
  void unget(char *line);
  uint64_t getUncPos()
  {
    return d_uncPos;
  }
  uint64_t uncompressedSize();
  void seek(uint64_t pos);
private:
  bool getChar(char* c);
  void skip(uint64_t toSkip);
  FILE* d_fp;
  
  struct ZState {
    ZState();
    ZState(const ZState& orig);
    ~ZState();
    ZState& operator=(const ZState& rhs);
    uint64_t fpos;
    z_stream s;
  } d_zs;
  int d_have;
  int d_datapos;

  char d_inbuffer[4096], d_outbuffer[32768];
  std::map<uint64_t, ZState> d_restarts;
  uint64_t d_uncPos;
  bool d_haveSeeked;
  std::string d_stash;
};

class BGZFWriter
{
public:
  BGZFWriter(const std::string& fname);
  ~BGZFWriter();
  uint64_t write(const char*, unsigned int len);

  void write32(uint32_t val);
  void writeBAMString(const std::string& str);
  void emitBlock(bool force=false);
private:

  void beginBlock();
  FILE* d_fp;
  z_stream d_s;
  std::string d_extra;
  gz_header d_gzheader;
  std::string d_block;
  uint32_t d_written;
  uint64_t d_blockstartpos;
};

void emitBGZF(FILE* fp, const std::string& block);
