#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include "zstuff.hh"
#include <stdexcept>
#include <iostream>
#include <boost/lexical_cast.hpp>

using namespace std;

ZLineReader::ZState::ZState()
{
  d_have=-1;
  d_datapos=-1;
  memset(&s, 0, sizeof(s));
  inflateInit2(&s, 31);
}
ZLineReader::ZState::~ZState()
{
  inflateEnd(&s);
}

ZLineReader::ZState::ZState(const ZLineReader::ZState& orig)
{
  *this=orig;
}


ZLineReader::ZState& ZLineReader::ZState::operator=(const ZLineReader::ZState& orig)
{
  memset(&s, 0, sizeof(s));
  fpos = orig.fpos;
  d_have = orig.d_have;
  d_datapos = orig.d_datapos;

  auto res=inflateCopy(&s, (z_stream*)&orig.s);
  if(res != Z_OK)  {
    throw std::runtime_error("Unable to copy Z state");
  }
  memcpy(d_inbuffer, orig.d_inbuffer, sizeof(d_inbuffer));
  memcpy(d_outbuffer, orig.d_outbuffer, sizeof(d_outbuffer));
  s.next_in = (Bytef*)(d_inbuffer + ((char*)orig.s.next_in - orig.d_inbuffer));
  s.avail_in = orig.s.avail_in;
  s.next_out=(Bytef*)d_outbuffer;
  s.avail_out=sizeof(d_outbuffer);
  return *this;
}


ZLineReader::ZLineReader(const std::string& fname)
{
  d_fp=fopen(fname.c_str(), "r");
  if(!d_fp)
    throw runtime_error("Unable to open '"+fname+"' for reading on ZLineReader"+ string(strerror(errno)));
  
  int ret = fread(d_zs.d_inbuffer, 1, sizeof(d_zs.d_inbuffer), d_fp);
  d_zs.s.avail_in=ret;
  d_zs.s.next_in = (Bytef*)d_zs.d_inbuffer;
  //  cerr<<"Got ret "<<ret<<endl;
  d_zs.s.next_out=(Bytef*)d_zs.d_outbuffer;
  d_zs.s.avail_out=sizeof(d_zs.d_outbuffer);
  d_zs.d_datapos=0;
  if(inflate(&d_zs.s, Z_NO_FLUSH) != Z_OK)
    throw runtime_error("Error inflating after open: "+string(d_zs.s.msg ? d_zs.s.msg : "no error message"));
  
  d_zs.d_have = d_zs.s.next_out - (Bytef*)d_zs.d_outbuffer;
  d_uncPos=0;
  d_haveSeeked=0;
}

bool ZLineReader::getChar(char* c)
{
  if(!d_zs.d_have) {
    //    cerr<<"Nothing available.."<<endl;
    d_zs.s.next_out=(Bytef*)d_zs.d_outbuffer;
    d_zs.s.avail_out=sizeof(d_zs.d_outbuffer);
    d_zs.d_datapos=0;
     
    if(d_zs.s.avail_in) {
      //      cerr<<"zlib has a bit more to chew on"<<endl;
      auto res=inflate(&d_zs.s, Z_NO_FLUSH);
      if(res!= Z_OK) {
        if(res == Z_STREAM_END) {
          if(d_zs.s.next_out == (Bytef*)d_zs.d_outbuffer)
            return false;
        }
        else
          throw runtime_error("Error inflating 1: "+ string(d_zs.s.msg ? d_zs.s.msg : "no error message"));
      }
      d_zs.d_have = d_zs.s.next_out - (Bytef*)d_zs.d_outbuffer;
    }
    if(!d_zs.d_have) {
      //      cerr<<"Still no output, getting more input.. "<<d_zs.s.avail_in<<endl;
      d_zs.s.next_in = (Bytef*)d_zs.d_inbuffer;
      d_zs.s.avail_in = fread(d_zs.s.next_in, 1, sizeof(d_zs.d_inbuffer) - d_zs.s.avail_in, d_fp);
      //      cerr<<"d_zs.s.avail_in: "<<d_zs.s.avail_in<<endl;
      if(!d_zs.s.avail_in)
        return false;
      auto res = inflate(&d_zs.s, Z_NO_FLUSH);
      if(res == Z_STREAM_END) {
        if(d_zs.s.next_out == (Bytef*)d_zs.d_outbuffer)
          return false;
      }
      else if(res != Z_OK)
        throw runtime_error("Error inflating 2: "+ string(d_zs.s.msg ? d_zs.s.msg : "no error message"));
      d_zs.d_have = d_zs.s.next_out - (Bytef*)d_zs.d_outbuffer;
    }
  }
  
  *c = d_zs.d_outbuffer[d_zs.d_datapos];
  //  cerr<<"Returning: '"<<*c<<"', d_zs.s.avail_out: "<<d_zs.d_have<<"\n";
  d_zs.d_datapos++;
  d_zs.d_have--;
  d_uncPos++;
  return true;
}

int ZLineReader::fgets(char* line, int num)
{
   if(!d_haveSeeked && (d_restarts.empty() || d_uncPos - d_restarts.rbegin()->first > 100000)) {
    d_zs.fpos = ftell(d_fp);
    d_restarts[d_uncPos]=d_zs;
  }

  char c;
  int i;
  for(i=0; i<num;++i) {
    if(!getChar(&c))
      break;
    *line++=c;
    if(c=='\n')
      break;
  }
  *line=0;

  return i;
}

bool ZLineReader::getLine(std::string* line)
{
  if(!d_haveSeeked && (d_restarts.empty() || d_uncPos - d_restarts.rbegin()->first > 100000)) {
    d_zs.fpos = ftell(d_fp);
    d_restarts[d_uncPos]=d_zs;
  }

  *line="";
  char buffer[80];
  char*p = buffer, *end=buffer+sizeof(buffer);
  char c;
  while(getChar(&c)) {
    if(c=='\n')
      break;
    *p++=c;
    if(p==end) {
      line->append(buffer, p);
      p=buffer;
    }
  }
  line->append(buffer, p);
  return !line->empty();
}

void ZLineReader::seek(uint64_t pos)
{
  d_haveSeeked=1;
  auto iter = d_restarts.lower_bound(pos);
  if(iter == d_restarts.end()) {
    throw runtime_error("Found nothing for pos = "+boost::lexical_cast<string>(pos));
  }
  if(iter != d_restarts.begin())
    --iter;
  //cerr<<"Want to seek to uncompressed pos: "<<pos<<", seeking to fpos: "<<iter->second.fpos;
  //cerr<<", giving us uncompressed pos "<<iter->first<<endl;
  fseek(d_fp, iter->second.fpos, SEEK_SET);
  d_zs = iter->second;
  d_uncPos = iter->first;
  char dontcare;
  //cerr<<"Now need to skip "<<pos - iter->first<<" bytes!"<<endl;
  for(unsigned int i = 0 ; i < pos - iter->first; ++i) {
    if(!getChar(&dontcare)) {
      throw runtime_error("Had EOF while seeking?!");
    }
  }
}

ZLineReader::~ZLineReader()
{
  inflateEnd(&d_zs.s);
  fclose(d_fp);
}

