#pragma once
#include <string>
#include <vector>

double getGCContent(const std::string& str);


template<typename T>
class Clusterer
{
public:
  explicit Clusterer(int limit) : d_limit(limit)
  {}

  void feed(const T& t)
  {
    if(d_clusters.empty() || t.pos - d_clusters.rbegin()->getEnd() > d_limit) {
      d_clusters.push_back(cluster());
    }
    d_clusters.rbegin()->d_members.push_back(t);
  }

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

    std::vector<T> d_members;
  };

  std::vector<cluster> d_clusters;
private:
  unsigned int d_limit;
};
