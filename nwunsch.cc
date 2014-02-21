#include <string.h>
#include <string>
#include <algorithm>
#include <iostream>

using std::min;
using std::swap;
using std::cout;
using std::endl;
using std::cerr;

void stringalign(const std::string& ain, const std::string& bin, double mispen, double gappen,
		 double skwpen, std::string& aout, std::string& bout, std::string& summary) {
  unsigned int i,j,k;
  double dn,rt,dg;
  std::string::size_type ia = ain.length(), ib = bin.length();
  aout.resize(ia+ib);
  bout.resize(ia+ib);
  summary.resize(ia+ib);
  double *cost[ia+1];
  for(unsigned int n=0; n < ia+1; ++n)
    cost[n]=new double[ib+1];
  cost[0][0] = 0.;
  for (i=1;i<=ia;i++) cost[i][0] = cost[i-1][0] + skwpen;
  for (i=1;i<=ib;i++) cost[0][i] = cost[0][i-1] + skwpen;
  for (i=1;i<=ia;i++) for (j=1;j<=ib;j++) {
      dn = cost[i-1][j] + ((j == ib)? skwpen : gappen);
      rt = cost[i][j-1] + ((i == ia)? skwpen : gappen);
      dg = cost[i-1][j-1] + ((ain[i-1] == bin[j-1])? -1. : mispen);
      cost[i][j] = min({dn,rt,dg});
    }
  i=ia; j=ib; k=0;
  while (i > 0 || j > 0) {
    dn = rt = dg = 9.99e99;
    if (i>0) dn = cost[i-1][j] + ((j==ib)? skwpen : gappen);
    if (j>0) rt = cost[i][j-1] + ((i==ia)? skwpen : gappen);
    if (i>0 && j>0) dg = cost[i-1][j-1] +
		      ((ain[i-1] == bin[j-1])? -1. : mispen);
    if (dg <= min(dn,rt)) {
      aout[k] = ain[i-1];
      bout[k] = bin[j-1];
      summary[k++] = ((ain[i-1] == bin[j-1])? '=' : '!');
      i--; j--;
    }
    else if (dn < rt) {
      aout[k] = ain[i-1];
      bout[k] = ' ';
      summary[k++] = ' ';		
      i--;
    }
    else {
      aout[k] = ' ';
      bout[k] = bin[j-1];
      summary[k++] = ' ';		
      j--;
    }
  }
  for (i=0;i<k/2;i++) {
    swap(aout[i],aout[k-1-i]);
    swap(bout[i],bout[k-1-i]);
    swap(summary[i],summary[k-1-i]);
  }
  aout.resize(k); bout.resize(k); summary.resize(k);
  for(unsigned int n=0; n < ia+1; ++n)
    delete[] cost[n];
}

int main(int argc, char**argv)
{
  if(argc!=3) {
    cerr<<"Syntax: nwunsch stringA stringB"<<endl;
    return EXIT_FAILURE;
  }

  std::string aout, bout, summary;
  stringalign(argv[1], argv[2], 1, 1, 0, aout, bout, summary);
  printf("A: %s\nB: %s\na: %s\nb: %s\nd: %s\n",
	 argv[1], argv[2],
	 aout.c_str(), bout.c_str(), summary.c_str());
}
