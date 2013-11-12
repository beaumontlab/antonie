#pragma once
#include <stdio.h>
#include <string>
#include <stdint.h>

void chomp(char* line);
char* sfgets(char* p, int num, FILE* fp);
void reverseNucleotides(std::string* nucleotides);
uint64_t filesize(const char* name);
bool stringfgets(FILE* fp, std::string* line);

struct VarMeanEstimator
{
  VarMeanEstimator() : N(0), xTot(0), x2Tot(0) {}
  void operator()(double val) 
  {
    ++N;
    xTot += val;
    x2Tot += val*val;
  }

  uint64_t N;
  double xTot;
  double x2Tot;
};

inline double mean(const VarMeanEstimator& vme)
{
  return vme.xTot/vme.N;
}

inline double variance(const VarMeanEstimator& vme) 
{
  return (vme.x2Tot - vme.xTot*vme.xTot/vme.N)/vme.N;
}
