#pragma once
#include <string>
#include <vector>
#include <functional>

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

