#pragma once
#include <string>
#include <vector>
#include <stdint.h>
#include <map>
#include <boost/flyweight.hpp>
#include <IntervalTree/IntervalTree.h>

//! A Gene annotation
struct GeneAnnotation
{
  boost::flyweight<std::string> chromosome;
  std::string tag;
  boost::flyweight<std::string> id;
  boost::flyweight<std::string> parent;
  boost::flyweight<std::string> type;
  std::string name;
  bool strand;
  uint64_t startPos;
  uint64_t stopPos;
  bool gene;
};

inline bool operator<(const GeneAnnotation&A, const GeneAnnotation& B)
{
  return A.startPos < B.startPos;
}

//! Provides GeneAnnotation objects as read from a GFF3 file
class GeneAnnotationReader
{
public:
  GeneAnnotationReader(const std::string& fname); //!< Parse GFF3 from fname
  std::vector<GeneAnnotation> lookup(std::string_view chromosome, uint64_t pos); //!< Get all annotations for pos
  std::vector<GeneAnnotation> lookup(std::string_view chromosome, uint64_t pos1, uint64_t pos2); //!< Get all annotations for pos

  std::vector<GeneAnnotation> getAll(std::string_view chromosome);
  uint64_t size() const
  {
    size_t ret{0};
    //    for(const auto& ga : d_gas)
      //      ret += ga.second.size();
    return ret;
  } //!< Number of annotations known

  std::vector<std::string> getChromosomes()
  {
    std::vector<std::string> ret;
    for(const auto& ga : d_gas)
      ret.push_back(ga.first);
    return ret;
  }

private:
  typedef IntervalTree<unsigned int, GeneAnnotation> gas_t;
  void parseGenBank(const std::string& fname);
  std::map<std::string, gas_t> d_gas;
};

std::vector<GeneAnnotation> parseGenBankString(const std::string& bank);
