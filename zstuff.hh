#pragma once
#include <string>
#include <zlib.h>
#include <stdio.h>
#include <map>
#include <stdexcept>

class ZLineReader
{
public:
  ZLineReader(const std::string& fname);
  ~ZLineReader();
  bool getLine(std::string* str);
  bool getChar(char* c);
  uint64_t getUncPos()
  {
    return d_uncPos;
  }
  void seek(uint64_t pos);
private:

  FILE* d_fp;
  
  struct ZState {
    ZState();
    ZState(const ZState& orig);
    uint64_t fpos;
    z_stream s;
    char d_inbuffer[4096], d_outbuffer[8192];
    int d_have;
    int d_datapos;

  } d_zs;
  std::map<uint64_t, ZState> d_restarts;
  uint64_t d_uncPos;
  bool d_haveSeeked;
};
