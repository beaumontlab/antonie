#include <iostream>
#include "misc.hh"
#include "fastq.hh"
#include <map>
using namespace std;

map<int,int> g_overlaps;

bool tryMerge(const FastQRead& one, const FastQRead& two, FastQRead* together)
{
  FastQRead inv(two);
  inv.reverse();
  
  if(inv.d_nucleotides.find(one.d_nucleotides.substr(0, 19)) != string::npos) {
    for(int overlap = one.d_nucleotides.length() ; overlap > 19; --overlap) {
      if(one.d_nucleotides.substr(0, overlap) == inv.d_nucleotides.substr(inv.d_nucleotides.length()-overlap)) {
	g_overlaps[overlap]++;
	//      cerr<<"Got overlap of "<<overlap<<endl;
	//      cerr<<one.d_nucleotides<<endl;
	//      cerr<<string(one.d_nucleotides.length()-overlap,' ')<<inv.d_nucleotides<<endl;
	together->d_nucleotides = inv.d_nucleotides;
	together->d_nucleotides = two.d_nucleotides.substr(overlap);

	//      cerr<<together->d_nucleotides<<endl;
	return true;
      }
    }
  }

  if(inv.d_nucleotides.find(one.d_nucleotides.substr(one.d_nucleotides.length()-19)) == string::npos)
    return false;

  for(int overlap = one.d_nucleotides.length() ; overlap > 19; --overlap) {
    if(one.d_nucleotides.substr(one.d_nucleotides.length()-overlap) == inv.d_nucleotides.substr(0, overlap)) {

      g_overlaps[overlap]++;
      //      cerr<<"Got overlap of "<<overlap<<endl;
      //      cerr<<one.d_nucleotides<<endl;
      //      cerr<<string(one.d_nucleotides.length()-overlap,' ')<<inv.d_nucleotides<<endl;
      together->d_nucleotides = one.d_nucleotides;
      together->d_nucleotides += inv.d_nucleotides.substr(overlap);
      //      cerr<<together->d_nucleotides<<endl;
      return true;
    }
  }
  
  return false;
}

bool findInRead(const FastQRead& fqr, const string& search, const string& rsearch)
{
  auto pos = fqr.d_nucleotides.find(search);
  if(pos != string::npos) {
    cout<<"F: "<<fqr.d_nucleotides.substr(pos)<<endl;
    return true;
  }
  else if((pos=fqr.d_nucleotides.find(rsearch)) != string::npos) {
    FastQRead inv = fqr;
    inv.reverse();
    pos = inv.d_nucleotides.find(search);
    cout<<"R: "<<inv.d_nucleotides.substr(pos)<<endl;
    return true;
  }
  return false;
}

int main(int argc, char**argv)
{
  string search(argv[1]);
  string rsearch(search);
  reverseNucleotides(&rsearch);

  StereoFASTQReader fqreader(argv[2], argv[3], 33);
  FastQRead fqr1, fqr2, merged;

  int mergedCount=0, total=0;
  while(fqreader.getReadPair(&fqr1, &fqr2)) {
    total++;
    if(tryMerge(fqr1, fqr2, &merged)) {
      mergedCount++;
      if(findInRead(merged, search, rsearch))
	cout<<"^ merged"<<endl;
    }
    else {
      findInRead(fqr1, search, rsearch);
      findInRead(fqr2, search, rsearch);
    }
  }
  cout<<"Got "<<total<<" pairs of which "<< mergedCount*100.0/total <<"% could be merged"<<endl;

  for(const auto& count : g_overlaps) {
    cout<<count.first<<'\t'<<count.second<<'\n';
  }
}


#if 0
#endif
