#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include "zstuff.hh"
#include <stdexcept>
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include "compat.hh"
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

using namespace std;

ZLineReader::ZState::ZState()
{
  fpos=0;
  memset(&s, 0, sizeof(s));
  inflateInit2(&s, 31);
}
ZLineReader::ZState::~ZState()
{
  inflateEnd(&s);
}

ZLineReader::ZState::ZState(const ZLineReader::ZState& orig)
{
  memset(&s, 0, sizeof(s));
  fpos = orig.fpos;

  auto res=inflateCopy(&s, (z_stream*)&orig.s);
  if(res != Z_OK)  {
    throw std::runtime_error("Unable to copy Z state");
  } 
}

ZLineReader::ZState& ZLineReader::ZState::operator=(const ZLineReader::ZState& orig)
{
  inflateEnd(&s);
  memset(&s, 0, sizeof(s));
  fpos = orig.fpos;

  auto res=inflateCopy(&s, (z_stream*)&orig.s);
  if(res != Z_OK)  {
    throw std::runtime_error("Unable to copy Z state");
  }
  return *this;
}

ZLineReader::ZLineReader(const std::string& fname)
{
  d_fp=fopen(fname.c_str(), "rb");
  if(!d_fp)
    throw runtime_error("Unable to open '"+fname+"' for reading on ZLineReader: "+ string(strerror(errno)));
  
  int ret = fread(d_inbuffer, 1, sizeof(d_inbuffer), d_fp);
  d_zs.s.avail_in=ret;
  d_zs.s.next_in = (Bytef*)d_inbuffer;
  //  cerr<<"Got ret "<<ret<<endl;
  d_zs.s.next_out=(Bytef*)d_outbuffer;
  d_zs.s.avail_out=sizeof(d_outbuffer);
  d_datapos=0;
  if(inflate(&d_zs.s, Z_NO_FLUSH) != Z_OK)
    throw runtime_error("Error inflating after open: "+string(d_zs.s.msg ? d_zs.s.msg : "no error message"));
  
  d_have = d_zs.s.next_out - (Bytef*)d_outbuffer;
  d_uncPos=0;
  d_haveSeeked=0;
}

bool ZLineReader::getChar(char* c)
{
  if(!d_have) {
    //    cerr<<"Nothing available.."<<endl;
    d_zs.s.next_out=(Bytef*)d_outbuffer;
    d_zs.s.avail_out=sizeof(d_outbuffer);
    d_datapos=0;
     
    if(d_zs.s.avail_in) {
      //      cerr<<"zlib has a bit more to chew on"<<endl;
      auto res=inflate(&d_zs.s, Z_NO_FLUSH);
      if(res!= Z_OK) {
        if(res == Z_STREAM_END) {
          if(d_zs.s.next_out == (Bytef*)d_outbuffer)
            return false;
        }
        else
          throw runtime_error("Error inflating 1: "+ string(d_zs.s.msg ? d_zs.s.msg : "no error message"));
      }
      d_have = d_zs.s.next_out - (Bytef*)d_outbuffer;
    }
    if(!d_have) {
      //      cerr<<"Still no output, getting more input.. "<<d_zs.s.avail_in<<endl;
      d_zs.s.next_in = (Bytef*)d_inbuffer;
      d_zs.s.avail_in = fread(d_zs.s.next_in, 1, sizeof(d_inbuffer) - d_zs.s.avail_in, d_fp);
      //      cerr<<"d_zs.s.avail_in: "<<d_zs.s.avail_in<<endl;
      if(!d_zs.s.avail_in)
        return false;
      auto res = inflate(&d_zs.s, Z_NO_FLUSH);
      if(res == Z_STREAM_END) {
        if(d_zs.s.next_out == (Bytef*)d_outbuffer)
          return false;
      }
      else if(res != Z_OK)
        throw runtime_error("Error inflating 2: "+ string(d_zs.s.msg ? d_zs.s.msg : "no error message"));
      d_have = d_zs.s.next_out - (Bytef*)d_outbuffer;
    }
  }
  
  if(c)
    *c = d_outbuffer[d_datapos];
  //  cerr<<"Returning: '"<<*c<<"', d_zs.s.avail_out: "<<d_have<<"\n";
  d_datapos++;
  d_have--;
  d_uncPos++;
  return true;
}

void ZLineReader::unget(char* line)
{
  d_stash=line;
}


char* ZLineReader::fgets(char* line, int num)
{
  if(!d_stash.empty()) {
    strncpy(line, d_stash.c_str(), num);
    d_stash.clear();
    return line;
  }
  if(!d_haveSeeked && (d_restarts.empty() || d_uncPos - d_restarts.rbegin()->first > 400000)) {
    d_zs.fpos = ftell(d_fp) - d_zs.s.avail_in;
    d_restarts[d_uncPos + d_have]=d_zs;
  }

  int i=0;
  char c;
  for(i=0; i<num;++i) {
    if(!getChar(&c))
      break;
    *line++=c;
    if(c=='\n')
      break;
  }
  *line=0;

  return i ? line : 0;
}

void ZLineReader::skip(uint64_t bytes)
{
  for(unsigned int i = 0 ; i < bytes; ++i) {
    if(!getChar(0)) {
      throw runtime_error("Had EOF while seeking?!");
    }
  }
}

uint64_t ZLineReader::uncompressedSize()
{
  throw runtime_error("Can't yet estimate size of compressed files, sorry");
}

void ZLineReader::seek(uint64_t pos)
{
  d_haveSeeked=1;
  
  auto iter = d_restarts.lower_bound(pos);
  if(iter == d_restarts.begin()) {
    throw runtime_error("Found nothing for pos = "+boost::lexical_cast<string>(pos));
  }
  if(iter != d_restarts.begin())
    --iter;
  //cerr<<"Want to seek to uncompressed pos: "<<pos<<", seeking to fpos: "<<iter->second.fpos;
  //cerr<<", giving us uncompressed pos "<<iter->first<<endl;


  if(pos > d_uncPos && (pos - d_uncPos) < (pos - iter->first)) {
    //    cerr<<"Skipping, "<< (pos - d_uncPos) << " < " << (pos - iter ->first)<<endl;
    skip(pos - d_uncPos);
    return;
  }

  fseek(d_fp, iter->second.fpos, SEEK_SET);
  d_zs = iter->second;
  d_uncPos = iter->first;

  d_have=0;
  d_datapos=0;
  d_zs.s.next_in=(Bytef*)d_inbuffer;
  d_zs.s.avail_in=0;
  d_zs.s.next_out=(Bytef*)d_outbuffer;
  d_zs.s.avail_out=sizeof(d_outbuffer);
  //cerr<<"Now need to skip "<<pos - iter->first<<" bytes!"<<endl;
  skip(pos - iter->first);

}

ZLineReader::~ZLineReader()
{
  inflateEnd(&d_zs.s);
  fclose(d_fp);
}

PlainLineReader::PlainLineReader(const std::string& fname)
{
  d_fp=fopen(fname.c_str(), "rb");
  if(!d_fp)
    throw runtime_error("Unable to open '"+fname+"' for reading on ZLineReader: "+ string(strerror(errno)));
  
}

void PlainLineReader::unget(char* line)
{
  d_stash=line;
}

char* PlainLineReader::fgets(char* line, int num)
{
  if(!d_stash.empty()) {
    strncpy(line, d_stash.c_str(), num);
    d_stash.clear();
    return line;
  }

  return ::fgets(line, num, d_fp);
}

void PlainLineReader::seek(uint64_t pos)
{
  if(fseek(d_fp, pos, SEEK_SET) < 0)
    throw runtime_error("Seeking in plain line reader: "+string(strerror(errno)));
}

uint64_t PlainLineReader::uncompressedSize()
{
  struct stat buf;
  if(fstat(fileno(d_fp), &buf) < 0) 
    throw runtime_error("Unable to determine size of file in PlainLineReader");
  
  return buf.st_size;
}


uint64_t PlainLineReader::getUncPos()
{
  return ftell(d_fp);
}

PlainLineReader::~PlainLineReader()
{
  fclose(d_fp);
}

unique_ptr<LineReader> LineReader::make(const std::string& fname)
{
  if(boost::ends_with(fname, ".gz"))
    return unique_ptr<LineReader>(new ZLineReader(fname));
  else
    return unique_ptr<LineReader>(new PlainLineReader(fname));
}


BGZFWriter::BGZFWriter(const std::string& fname)
{
  if(fname.empty()) {
    d_fp=0;
    return;
  }
  d_fp=fopen(fname.c_str(), "w");
  if(!d_fp)
    throw runtime_error("Unable to open '"+fname+"' for BGZFWriter: "+ string(strerror(errno)));
  beginBlock();
}

void BGZFWriter::beginBlock()
{
  d_blockstartpos = ftell(d_fp);
  d_block.clear();
  memset(&d_gzheader, 0, sizeof(d_gzheader));
  d_gzheader.time=0;
  d_gzheader.os=0xff;
  d_extra.assign(1, 66);
  d_extra.append(1, 67);
  d_extra.append(1, 2);
  d_extra.append(1, (char)0);

  d_extra.append(2, 0); // our size, we spoof it in afterwards
  
  d_gzheader.extra=(Bytef*)d_extra.c_str();
  d_gzheader.extra_len=d_extra.length();

  memset(&d_s, 0, sizeof(d_s));
  auto res = deflateInit2(&d_s, Z_DEFAULT_COMPRESSION, Z_DEFLATED, 16+15, 8, Z_DEFAULT_STRATEGY);
  if(res != Z_OK) {
    throw runtime_error("Unable to initialize compression");
  }
  
  if(deflateSetHeader(&d_s, &d_gzheader) != Z_OK)
    throw runtime_error("Unable to initialize compression header");

  d_written=0;
}

uint64_t BGZFWriter::write(const char*c, unsigned int len)
{
  uint64_t pos = (d_blockstartpos<<16) | (d_written);
  d_s.next_in = (Bytef*) c;
  d_s.avail_in = len;

  char *buffer = new char[1024+len*2]; // XXX windows
  do {
    d_s.next_out = (Bytef*) buffer;
    d_s.avail_out=sizeof(buffer);

    //cerr<<d_s.avail_out<<endl;
    //cerr<<"avail_in: "<<d_s.avail_in<<", flush: "<<flush<<endl;
    auto res = deflate(&d_s, 0);
    //cerr<<"avail_in post: "<<d_s.avail_in<<", res="<<res<<endl;
    if(res != Z_OK)
      throw runtime_error("Unable to deflate data: "+string(d_s.msg ? d_s.msg : "no error"));
    //cerr<<"Got "<<(d_s.next_out - (Bytef*)buffer)<<" bytes (end="<<(res==Z_STREAM_END)<<")"<<endl;
    d_block.append(buffer, d_s.next_out - (Bytef*)buffer);
  } while(d_s.avail_in);
  delete[] buffer;
  d_written+=len;
  if(d_written > 65000)
    emitBlock();
  return pos;
}

void BGZFWriter::emitBlock(bool force)
{
  if(!d_written && !force)
    return;
  char buffer[65535];
  d_s.next_out = (Bytef*) buffer;
  d_s.avail_out=sizeof(buffer);
  d_s.avail_in=0;
  auto res = deflate(&d_s, Z_FINISH);
  if(res != Z_STREAM_END) {
    throw runtime_error("Shit");
  }
  auto got = d_s.next_out - (Bytef*) buffer;
  d_block.append(buffer, got);
  uint16_t len = d_block.size()-1;
  d_block.replace(16, 2, (char*)&len, 2);
  fwrite(d_block.c_str(), 1, d_block.size(), d_fp);
  //  cerr<<"Writing block of "<<d_block.size()<<" ("<<d_written<<")"<<endl;
  d_block.clear();
  d_extra.clear();
  d_written=0;
  deflateEnd(&d_s);
  beginBlock();
}

void BGZFWriter::write32(uint32_t val)
{
  uint32_t nval = htonl(val); 
  write((char*)&nval, 4);
}

void BGZFWriter::writeBAMString(const std::string& str)
{
  write32(str.length());
  write(str.c_str(), str.length());
}

BGZFWriter::~BGZFWriter()
{
  if(d_fp==0)
    return;
  emitBlock();
  emitBlock(true);
  deflateEnd(&d_s);
  fclose(d_fp);
}
