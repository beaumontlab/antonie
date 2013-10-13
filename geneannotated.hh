#pragma once
#include <string>
#include <vector>

struct GeneAnnotation
{
  std::string tag;
  std::string name;
  uint64_t startPos;
  uint64_t stopPos;
};

inline bool operator<(const GeneAnnotation&A, const GeneAnnotation& B)
{
  return A.startPos < B.startPos;
}

class GeneAnnotationReader
{
public:
  GeneAnnotationReader(const std::string& fname);
  std::vector<GeneAnnotation> lookup(uint64_t pos);
  uint64_t size() const { return d_gas.size(); }
private:
  typedef std::vector<GeneAnnotation> gas_t;
  gas_t d_gas;
};

