#include "misc.hh"
#include <string.h>
#include <stdexcept>


char* sfgets(char* p, int num, FILE* fp)
{
  char *ret = fgets(p, num, fp);
  if(!ret) 
    throw std::runtime_error("Unexpected EOF");
  return ret;
}

void chomp(char* line)
{
  char *p;
  p = strchr(line, '\r');
  if(p)*p=0;
  p = strchr(line, '\n');
  if(p)*p=0;
}
