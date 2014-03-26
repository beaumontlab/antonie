// read all input files, output for each line of text in which of the input files it was found
// public domain code by bert.hubert@netherlabs.nl
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <boost/dynamic_bitset.hpp>
#include <boost/algorithm/string.hpp>

using namespace std;

// this allows us to make a Case Insensitive container
struct CIStringCompare: public std::binary_function<string, string, bool> 
{
  bool operator()(const string& a, const string& b) const
  {
    if(std::all_of(a.begin(), a.end(), ::isdigit) && 
       std::all_of(b.begin(), b.end(), ::isdigit))
      return atoi(a.c_str()) < atoi(b.c_str());
    return strcasecmp(a.c_str(), b.c_str()) < 0;
  }
};

int main(int argc, char**argv)
{
  typedef map<string, boost::dynamic_bitset<>, CIStringCompare> presence_t;
  presence_t presence;

  string line;
  cout << '\t';
  for(int n = 1; n < argc; ++n) {
    cout << argv[n] << '\t';
    ifstream ifs(argv[n]);
    if(!ifs) {
      cerr<<"Unable to open '"<<argv[n]<<"' for reading\n"<<endl;
      exit(EXIT_FAILURE);
    }

    while(getline(ifs, line)) {
      boost::trim(line);
      if(line.empty())
	continue;
      presence_t::iterator iter = presence.find(line);
      if(iter == presence.end()) { // not present, do a very efficient 'insert & get location'
	iter = presence.insert(make_pair(line, boost::dynamic_bitset<>(argc-1))).first; 
      }
      iter->second[n-1]=1;
    }
  }
  cout << '\n';

  // this is where we store the reverse map, 'presence groups', so which lines where present in file1, but not file2 etc
  typedef map<boost::dynamic_bitset<>, vector<string> > revpresence_t;
  revpresence_t revpresence;

  for(const auto& val : presence) {
    revpresence[val.second].push_back(val.first);
    cout << val.first << '\t';
    for (boost::dynamic_bitset<>::size_type i = 0; i < val.second.size(); ++i) {
      cout << val.second[i] << '\t';
    }
    cout << endl;
  }

  cout << "\nPer group output\t\n";
  for(const auto& val : revpresence) {
    cout<<"\nGroup: \t";
    for (boost::dynamic_bitset<>::size_type i = 0; i < val.first.size(); ++i) {
      cout << val.first[i] << '\t';
    }
    cout << endl << "       \t";
    for (boost::dynamic_bitset<>::size_type i = 0; i < val.first.size(); ++i) {
      if(val.first[i])
	cout<<argv[i+1];
      cout<<'\t';
    }
    cout<<endl;

    for(const auto& entry : val.second) {
      cout << entry << "\t\n";
    }
  }
}
