#pragma once
#include <string>
#include <vector>
#include <functional>
#include <stdlib.h>
#include <sstream> 
#include "antonie.hh"

extern const char* g_gitHash;

//! convert a Sanger Q-score into an error probability. Uses a cache to be fast.
double qToErr(unsigned int i);

//! returns GC fraction of nucleotides in str
double getGCContent(const std::string& str);


//! Generic class to cluster objects that are 'close by'
template<typename T>
class Clusterer
{
public:
  //! Cluster objects that are less than 'limit' apart together
  explicit Clusterer(int limit) : d_limit(limit)
  {}

  //! Feed an object
  void feed(const T& t)
  {
    if(d_clusters.empty() || t.pos - d_clusters.rbegin()->getEnd() > d_limit) {
      d_clusters.push_back(cluster());
    }
    d_clusters.rbegin()->d_members.push_back(t);
  }

  //! Represents a cluster
  struct cluster
  {
    int getBegin()
    {
      return d_members.begin()->pos;
    }
    int getEnd()
    {
      return d_members.rbegin()->pos;
    }
    int getMiddle()
    {
      return (getBegin()+getEnd())/2;
    }

    std::vector<T> d_members; //!< members of this cluster
  };

  //! The clusters we made for you
  std::vector<cluster> d_clusters;
private:
  unsigned int d_limit;
};

typedef std::function<void(void)> acgt_t;
inline void acgtDo(char c, acgt_t aDo, acgt_t cDo, acgt_t gDo, acgt_t tDo)
{
  switch(c) {
  case 'A':
    aDo();
    break;
  case 'C':
    cDo();
    break;
  case 'G':
    gDo();
    break;
  case 'T':
    tDo();
    break;
  }
}

typedef std::function<void(void)> acgt_t;
inline void acgtxDo(char c, acgt_t aDo, acgt_t cDo, acgt_t gDo, acgt_t tDo, acgt_t xDo)
{
  switch(c) {
  case 'A':
    aDo();
    break;
  case 'C':
    cDo();
    break;
  case 'G':
    gDo();
    break;
  case 'T':
    tDo();
    break;
  case 'X':
    xDo();
    break;
  }
}


//! Little utility to pick a random element from a container
template<typename T>
const typename T::value_type& pickRandom(const T& t)
{
  return t[rand() % t.size()];
}

/** returns v as a string in JSON format, where v is a vector of values, and we return it as an array of [offset,value] pairs. v can be transformed inline  ysing yAdjust and xAdjust
    \param v Vector of values
    \param name name of JSON object
    \param yAdjust function (or lambda) that transforms the values in v
    \param xAdjust function (or lambda) that generates the offsets in our return vector. Gets passed this offset, returns a double
*/
template<typename T>
std::string jsonVector(const std::vector<T>& v, const std::string& name, 
		  std::function<double(dnapos_t)> yAdjust = [](dnapos_t d){return 1.0*d;},
		  std::function<double(int)> xAdjust = [](int i){return 1.0*i;})
{
  std::ostringstream ret;
  ret << "var "<<name<<"=[";
  for(auto iter = v.begin(); iter != v.end(); ++iter) {
    if(iter != v.begin())
      ret<<',';
    ret << '[' << xAdjust(iter - v.begin()) <<','<< yAdjust(*iter)<<']';
  }
  ret <<"];\n";
  return ret.str();
}


template<typename T>
std::string jsonVectorD(const std::vector<T>& v, 
		  std::function<double(double)> yAdjust = [](double d){return d;},
		  std::function<double(int)> xAdjust = [](int i){return 1.0*i;})
{
  std::ostringstream ret;
  ret<<"[";
  for(auto iter = v.begin(); iter != v.end(); ++iter) {
    if(iter != v.begin())
      ret<<',';
    ret << '[' << std::fixed<< xAdjust(iter - v.begin()) <<','<< std::scientific << yAdjust(*iter)<<']';
  }
  ret <<"]";
  return ret.str();
}
