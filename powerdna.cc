#include <stdio.h>
#include <iostream>
#include <stdint.h>
#include <map>
#include <vector>
#include <string>
#include <string.h>
#include <stdexcept>
#include <boost/foreach.hpp>
#include <boost/progress.hpp>
#include <algorithm>
#include <numeric>
#include <unistd.h>
#include <errno.h>
#include <math.h>
#include <boost/format.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/density.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/tuple/tuple.hpp>
#include "geneannotated.hh"
#include "misc.hh"
#include "fastq.hh"

#include <mba/diff.h>
#include <mba/msgno.h>


using namespace boost;
using namespace boost::accumulators;

extern "C" {
#include "hash.h"
}
using namespace std;

struct FASTQMapping
{
  uint64_t pos;
  bool reverse;
  int indel; // 0 = nothing, >0 means WE have an insert versus reference at pos
             // <0 means WE have a delete versus reference at pos
};

struct GenomeLocusMapping
{
  GenomeLocusMapping() : coverage(0) {}
  vector<FASTQMapping> d_fastqs;
  unsigned int coverage;
};
 
class ReferenceGenome
{
public:
  ReferenceGenome(const string& fname);
  uint64_t size() const {
    return d_genome.size();
  }

  vector<uint64_t> getReadPositions(const std::string& nucleotides)
  {
    vector<uint64_t> ret;
    if(nucleotides.length() != d_indexlength)
      throw runtime_error("Attempting to find fragment of length we've not indexed for");
      
    uint32_t hashval = hash(nucleotides.c_str(), nucleotides.length(), 0);
    HashPos hp(hashval, 0);
    pair<index_t::const_iterator, index_t::const_iterator> range = equal_range(d_index.begin(), d_index.end(), hp);
    if(range.first == range.second)
      return ret;

    for(;range.first != range.second; range.first++) {
      if(!memcmp(d_genome.c_str() + range.first->d_pos, nucleotides.c_str(), nucleotides.length())) {
	ret.push_back(range.first->d_pos);
      }
    }
    return ret;
  }
  
  uint64_t getReadPosBoth(FastQRead* fq) // tries both
  {
    vector<uint64_t> positions;
    string nucleotides;
    for(int tries = 0; tries < 2; ++tries) {
      positions = getReadPositions(fq->d_nucleotides);

      if(!positions.empty()) {
	if(positions.size() > 1)
	  random_shuffle(positions.begin(), positions.end());

	cover(positions[0], fq->d_nucleotides.size(), fq->d_quality);
	return positions[0];
      }
      fq->reverse();
    }
    return string::npos;
  }

  void cover(uint64_t pos, unsigned int length, const std::string& quality) 
  {
    const char* p = quality.c_str();
    for(unsigned int i = 0; i < length; ++i) {
      if(p[i]-33 > 30)
	d_mapping[pos+i].coverage++;
    }
  }

  void cover(uint64_t pos, char quality) 
  {
    if(quality-33 > 30)
      d_mapping[pos].coverage++;
  }

  void mapFastQ(uint64_t pos, const FastQRead& fqfrag, int indel=0)
  {
    FASTQMapping fqm;
    fqm.pos=fqfrag.position;
    fqm.reverse = fqfrag.reversed;
    fqm.indel = indel;
    d_mapping[pos].d_fastqs.push_back(fqm);
  }

  string snippet(uint64_t start, uint64_t stop) const { 
    return d_genome.substr(start, stop-start);
  }
  void printCoverage();
  void index(int length);
  void printFastQs(uint64_t pos, FASTQReader& fastq);
  vector<GenomeLocusMapping> d_mapping;
private:
  unsigned int d_indexlength;
  string d_genome;

  struct HashPos {
    HashPos(uint32_t hash_, uint64_t pos) : d_hash(hash_), d_pos(pos)
    {}
    HashPos(){}
    uint32_t d_hash;
    uint64_t d_pos;
    
    bool operator<(const HashPos& rhs) const 
    {
      return d_hash < rhs.d_hash;
    }
  };

  typedef vector<HashPos> index_t;
  index_t d_index;
};

ReferenceGenome::ReferenceGenome(const string& fname)
{
  FILE* fp = fopen(fname.c_str(), "r");
  if(!fp)
    throw runtime_error("Unable to open reference genome file '"+fname+"'");
  d_genome.reserve(filesize(fname.c_str()));  // slight overestimate which is great
  
  char line[256]="";

  sfgets(line, sizeof(line), fp);
  chomp(line);

  if(line[0] != '>') 
    throw runtime_error("Input not FASTA");
  cerr<<"Reading FASTA reference genome of '"<<line+1<<"'\n";

  d_genome="*"; // this gets all our offsets ""right""
  while(fgets(line, sizeof(line), fp)) {
    chomp(line);
    d_genome.append(line);
  }
  
  d_mapping.resize(d_genome.size());
  d_indexlength=0;
}

void ReferenceGenome::index(int length)
{
  d_index.clear();
  d_index.reserve(d_genome.length());
  d_indexlength=length;
  cerr<<"Indexing "<<d_genome.length()<<" nucleotides for length "<<length<<"..";

  //  boost::progress_display show_progress( d_genome.length() - d_indexlength, cerr);

  for(string::size_type pos = 0 ; pos < d_genome.length() - length; ++pos) {
    uint32_t hashval = hash(d_genome.c_str() + pos, d_indexlength, 0);
    d_index.push_back(HashPos(hashval, pos));
    // ++show_progress;
  }

  cerr<<" done\nSorting hashes..";
  sort(d_index.begin(), d_index.end());
  cerr<<" done"<<endl;
    
  //  cerr<<"Average hash fill: "<<1.0*d_genome.length()/d_index.size()<<endl;
}

void ReferenceGenome::printFastQs(uint64_t pos, FASTQReader& fastq) 
{
  if(pos < 150)
    pos = 0;
  else
    pos -= 150;
  
  string reference=snippet(pos, pos+300);
  unsigned int insertPos=0;
  for(unsigned int i = 0 ; i < 151; ++i) {
    if(i==75)
      cout << reference << endl;
    string spacer(i, ' ');
    BOOST_FOREACH(FASTQMapping& fqm, d_mapping[pos+i].d_fastqs) {
      fastq.seek(fqm.pos);
      FastQRead fqr;
      fastq.getRead(&fqr);
      if(fqm.reverse)
	fqr.reverse();

      if(fqm.indel > 0 && !insertPos) { // our read has an insert at this position
	reference.insert(i+fqm.indel, 1, '_');
	insertPos=i+fqm.indel;
	//	fqr.d_nucleotides.erase(fqm.indel, 1); // this makes things align again
	//fqr.d_quality.erase(fqm.indel, 1); 
      } else if(fqm.indel < 0) {      // our read has an erase at this position
	fqr.d_nucleotides.insert(-fqm.indel, 1, 'X');
	fqr.d_quality.insert(-fqm.indel, 1, 'X');
      }
      
      if(fqm.indel <= 0 && insertPos && i > insertPos) {
	fqr.d_nucleotides.insert(0, 1, '<');
	fqr.d_quality.insert(0, 1, 'X');
      }
      cout << spacer;
      int offset=0;
      for(unsigned int j = 0 ; j < fqr.d_nucleotides.size(); ++j) {
	if(reference[i+j]=='_' && !fqm.indel) {
	  cout<<'_';
	  offset=1;
	}
	if(reference[i+j+offset]==fqr.d_nucleotides[j])
	  cout<<'.';
	else if(fqr.d_quality[j] > '@') 
	  cout << fqr.d_nucleotides[j];
	else if(fqr.d_quality[j] < '7') 
	  cout << ' ';
	else
	  cout<< (char)tolower(fqr.d_nucleotides[j]);
      }
      cout<<"                 "<<(fqm.reverse ? 'R' : ' ');
      cout<<endl;
    }
  }
}

struct Unmatched
{
  string left, unmatched, right;
  uint64_t pos;
};

vector<Unmatched> g_unm;

void ReferenceGenome::printCoverage()
{
  uint64_t totCoverage=0, noCoverages=0;
  unsigned int cov;

  bool wasNul=true;
  string::size_type prevNulpos=0;

  vector<unsigned int> covhisto;
  covhisto.resize(65535);

  for(string::size_type pos = 0; pos < d_mapping.size(); ++pos) {
    cov = d_mapping[pos].coverage;
    bool noCov = cov < 1;
    covhisto[cov]++;
    totCoverage += cov;

    if(noCov) {
      noCoverages++;

    }
    
    if(!noCov && wasNul) {
      if(prevNulpos > 40 && pos + 40 < d_genome.length()) {
	Unmatched unm;
	unm.left = d_genome.substr(prevNulpos-40, 40);
	unm.right = d_genome.substr(pos, 40);
	unm.unmatched = d_genome.substr(prevNulpos, pos-prevNulpos);
	unm.pos = prevNulpos;
	g_unm.push_back(unm);
      }
      wasNul=false;
    }
    else if(noCov && !wasNul) {
      wasNul=true;
      prevNulpos = pos;
    }
  }

  cerr << (boost::format("Average depth: %|40t|    %10.2f\n") % (1.0*totCoverage/d_mapping.size())).str();
  cerr << (boost::format("Uncovered nucleotides: %|40t| %10d (%.2f%%)\n") % noCoverages % (noCoverages*100.0/d_mapping.size())).str();

  uint64_t total = std::accumulate(covhisto.begin(), covhisto.end(), 0);

  cout<<"var histo=[";
  for(unsigned int i = 0; i < 250; i++ ) 
  {
    if(i)
      cout<<',';
    std::cout << '['<< i << ',' << 1.0*covhisto[i]/total << ']';
  }
  cout<<"];\n";
}


struct LociStats
{
  vector<boost::tuple<char,char,char>> samples;
};

typedef map<uint64_t, LociStats> locimap_t;
locimap_t locimap;

// 0 if nothing interesting, positive if our read has insert at that position, negative if we have a delete at that position
int MBADiff(uint64_t pos, const FastQRead& fqr, const string& reference)
{
  string::size_type n, m, d;
  int sn, i;
  struct varray *ses = varray_new(sizeof(struct diff_edit), NULL);
  
  n = reference.length();
  m = fqr.d_nucleotides.length();
  if ((d = diff(reference.c_str(), 0, n, fqr.d_nucleotides.c_str(), 0, m, NULL, NULL, NULL, 0, ses, &sn, NULL)) == -1) {
    MMNO(errno);
    printf("Error\n");
    return EXIT_FAILURE;
  }
  int ret=0;
  if(sn == 4 && d == 2) {
    struct diff_edit *match1 = (struct diff_edit*)varray_get(ses, 0),
      *change1=(struct diff_edit*)varray_get(ses, 1), 
      *match2=(struct diff_edit*)varray_get(ses, 2), 
      *change2=(struct diff_edit*)varray_get(ses, 3);
    
    if(match1->op == DIFF_MATCH && match2->op==DIFF_MATCH) {
      if(change1->op == DIFF_DELETE && change2->op==DIFF_INSERT) {
	cout << "Have delete of "<<change1->len<<" in our read at " << pos+change1->off <<endl;
	ret=-change1->off;
      }
      else if(change1->op == DIFF_INSERT && change2->op==DIFF_DELETE) {
	cout<<"Have insert of "<<change1->len<<" in our read at "<<pos+change1->off<<endl;
	ret=change1->off;
      }
    }
  }
  if(sn > 6)
    return 0;

  printf("pos %ld, d=%ld sn=%d\nUS:  %s\nREF: %s\n", pos, d, 
	 sn, fqr.d_nucleotides.c_str(), reference.c_str());

  // should have 4 components, a MATCH, an INSERT/DELETE followed by a MATCH, followed by a DELETE/INSERT (inverse)

  for (i = 0; i < sn; i++) {
    struct diff_edit *e = (struct diff_edit*)varray_get(ses, i);

    switch (e->op) {
    case DIFF_MATCH:
      printf("MAT: ");
      fwrite(fqr.d_nucleotides.c_str() + e->off, 1, e->len, stdout);
      break;
    case DIFF_INSERT:
      printf("INS: ");
      fwrite(reference.c_str() + e->off, 1, e->len, stdout);
      break;
    case DIFF_DELETE:
      printf("DEL: ");
      fwrite(fqr.d_nucleotides.c_str() + e->off, 1, e->len, stdout);
      break;
    }
    printf("\n");
  }
  varray_del(ses);
  return ret;
}

string DNADiff(ReferenceGenome& rg, uint64_t pos, FastQRead& fqfrag)
{
  string diff;
  diff.reserve(fqfrag.d_nucleotides.length());

  string reference = rg.snippet(pos, pos + fqfrag.d_nucleotides.length());

  int diffcount=0;
  for(string::size_type i = 0; i < fqfrag.d_nucleotides.size() && i < reference.size();++i) {
    if(fqfrag.d_nucleotides[i] != reference[i] && fqfrag.d_quality[i]>'@') 
      diffcount++;
  }

  if(diffcount < 5) 
    rg.mapFastQ(pos, fqfrag);
  else {
    int res=MBADiff(pos, fqfrag, reference);
    if(res) {
      rg.mapFastQ(pos, fqfrag, res);
      diffcount=1;
      if(res > 0) { // our read has an insert at this position
	fqfrag.d_nucleotides.erase(res, 1); // this makes things align again
	fqfrag.d_quality.erase(res, 1); 
      } else {      // our read has an erase at this position
	fqfrag.d_nucleotides.insert(-res, 1, 'X');
	fqfrag.d_quality.insert(-res, 1, 'X');
      }
    }
  }
  //  cout<<"US:  "<<fqfrag.d_nucleotides<<endl<<"DIF: ";
  for(string::size_type i = 0; i < fqfrag.d_nucleotides.size() && i < reference.size();++i) {
    if(fqfrag.d_nucleotides[i] != reference[i]) {
      diff.append(1, fqfrag.d_quality[i] > '@' ? '!' : '^');
      if(fqfrag.d_quality[i]>'@' && diffcount < 5 && ( (!fqfrag.reversed && i > 14) || (fqfrag.reversed && i < (150-14) ) ) ) 
	locimap[pos+i].samples.push_back(boost::make_tuple(fqfrag.d_nucleotides[i], fqfrag.d_quality[i], fqfrag.reversed ^ (i>75))); // head or tail
    }
    else {
      diff.append(1, ' ');
      rg.cover(pos+i,fqfrag.d_quality[i]);
    }
  }

  //  cout<<diff<<endl<<"REF: "<<b<<endl;
  //  cout<<"QUA: "<<fqfrag.d_quality<<endl;
  return diff;
}


void printUnmatched(ReferenceGenome& rg, const string& name)
{
  unsigned int unmcount=0;
  BOOST_FOREACH(Unmatched& unm, g_unm) {
    printf("%s[%d]=[", name.c_str(), unmcount++);
    for(uint64_t pos = unm.pos - 500; pos < unm.pos + 500; ++pos) {
      if(pos != unm.pos - 500) 
	printf(", ");
      printf("[%ld, %d]", pos, rg.d_mapping[pos].coverage);
    }
    printf("];\n");
  }
  cout<<endl;
}

unsigned int variabilityCount(const ReferenceGenome& rg, uint64_t position, const LociStats& lc, double* fraction)
{
  vector<int> counts(256);
  counts[rg.snippet(position, position+1)[0]]+=rg.d_mapping[position].coverage;
  
  int forwardCount=0;

  for(vector<boost::tuple<char,char,char> >::const_iterator j = lc.samples.begin(); j!= lc.samples.end(); ++j) {
    counts[j->get<0>()]++;
    if(j->get<2>())
      forwardCount++;
  }
  sort(counts.begin(), counts.end());
  unsigned int nonDom=0;
  for(unsigned int i=0; i < 255; ++i) {
    nonDom+=counts[i];
  }
  
  if(nonDom + counts[255] < 20)
    return 0;

  *fraction = 1.0*forwardCount / (1.0*lc.samples.size());
  if(*fraction < 0.05 || *fraction > 0.9)
    return 0;

  return 100*nonDom/counts[255];
}

int main(int argc, char** argv)
{
  if(argc!=3) {
    fprintf(stderr, "Syntax: powerdna fastqfile fastafile\n");
    exit(EXIT_FAILURE);
  }
  GeneAnnotationReader gar("NC_012660.csv");
  cerr<<"Done reading "<<gar.size()<<" annotations"<<endl;
  
  FASTQReader fastq(argv[1]);

  unsigned int bytes=0;
  FastQRead fqfrag;
  bytes=fastq.getRead(&fqfrag); // get a fragment to index based on its size
  
  ReferenceGenome rg(argv[2]);

  rg.index(fqfrag.d_nucleotides.size());

  cerr<<"Loading positive control filter genome(s)"<<endl;
  ReferenceGenome phix("phi-x174.fasta");
  phix.index(fqfrag.d_nucleotides.size());

  uint64_t pos;

  uint64_t withAny=0, found=0, notFound=0, total=0, qualityExcluded=0, fuzzyFound=0, phixFound=0;

  cerr<<"Performing exact matches of reads to reference genome"<<endl;
  boost::progress_display show_progress( filesize(argv[1]), cerr);

  accumulator_set<double, stats<tag::mean, tag::median > > acc;

  vector<uint64_t> unfoundReads;
  do { 
    show_progress += bytes;
    total++;
    BOOST_FOREACH(char c, fqfrag.d_quality) {
      double i = c-33;
      acc(i);
    }

    if(fqfrag.d_nucleotides.find('N') != string::npos) {
      withAny++;
      continue;
    }
    pos = rg.getReadPosBoth(&fqfrag);
    if(pos == string::npos ) {
      if(phix.getReadPosBoth(&fqfrag)!=string::npos) {
	phixFound++;
      }
      else {
	unfoundReads.push_back(fqfrag.position);
	notFound++;
      }
    }
    else {
      rg.mapFastQ(pos, fqfrag);
      found++;
    }
  } while((bytes=fastq.getRead(&fqfrag)));

  cerr<< (boost::format("Total fragments: %|40t| %10d") % total).str() <<endl;
  cerr<< (boost::format("Excluded control fragments: %|40t|-%10d") % phixFound).str() <<endl;
  cerr<< (boost::format("Quality excluded: %|40t|-%10d") % qualityExcluded).str() <<endl;
  cerr<< (boost::format("Ignored reads with N: %|40t|-%10d") % withAny).str()<<endl;
  cerr<< (boost::format("Full matches: %|40t|-%10d (%.02f%%)\n") % found % (100.0*found/total)).str();
  cerr<< (boost::format("Not found: %|40t|=%10d (%.02f%%)\n") % notFound % (notFound*100.0/total)).str();

  cerr << (boost::format("Mean Q: %|40t|    %10.2f\n") % mean(acc)).str();
  cerr << (boost::format("Median Q: %|40t|    %10.2f\n") % median(acc)).str();

  rg.printCoverage();
  rg.printFastQs(45881, fastq);  
  
  printUnmatched(rg,"perfect");
  cout<<"---------"<<endl;

  int keylen=11;
  rg.index(keylen);

  string left, middle, right;
  vector<uint64_t> lpositions, mpositions, rpositions;
  vector<uint64_t> stillUnfound;
  typedef pair<uint64_t, char> tpos;
  cerr<<"Performing sliding window partial matches"<<endl;
  boost::progress_display fuzzyProgress(unfoundReads.size(), cerr);
  
  BOOST_FOREACH(uint64_t pos, unfoundReads) {
    fastq.seek(pos);
    fastq.getRead(&fqfrag);

    for(unsigned int attempts=0; attempts < 12; ++attempts) {
      for(int tries = 0; tries < 2; ++tries) {
	if(tries)
	  fqfrag.reverse();
	left=fqfrag.d_nucleotides.substr(attempts*3, keylen);	middle=fqfrag.d_nucleotides.substr(50+attempts*3, keylen); right=fqfrag.d_nucleotides.substr(100+attempts*3, keylen);
      	lpositions=rg.getReadPositions(left); 
	if(lpositions.empty())
	  continue;
	mpositions=rg.getReadPositions(middle); 
	if(mpositions.empty())
	  continue;
	rpositions=rg.getReadPositions(right);

	if(lpositions.size() + mpositions.size() + rpositions.size() < 3)
	  continue;
       	vector<tpos> together;
	BOOST_FOREACH(uint64_t fpos, lpositions) { together.push_back(make_pair(fpos, 'L')); }	
	BOOST_FOREACH(uint64_t fpos, mpositions) { together.push_back(make_pair(fpos, 'M')); }
	BOOST_FOREACH(uint64_t fpos, rpositions) { together.push_back(make_pair(fpos, 'R')); }
	
	sort(together.begin(), together.end());

	for(unsigned int i = 0; i < together.size() - 2; ++i) {
	  if(together[i].second=='L' && together[i+1].second=='M' && together[i+2].second=='R' && 
	     together[i+1].first - together[i].first < 60 && together[i+2].first - together[i+1].first < 60) {
	    uint64_t lpos, mpos, rpos;
	    lpos=together[i].first; mpos=together[i+1].first; rpos=together[i+2].first;
	    /*	    cout<<lpos<<together[i].second<<" - " << mpos<<together[i+1].second<<
	      " - " << rpos<<together[i+2].second<< " (" << 
	      mpos - lpos<< " - " << rpos - mpos <<"), reversed: "<<tries<<endl;
	    */
	    if(lpos < 3*attempts)
	      continue;
	    
	    uint64_t matchOffset=lpos - 3*attempts;
	    DNADiff(rg, matchOffset, fqfrag);
	    fuzzyFound++;
	    goto foundIt;
	    break;
	  }
	}
      }
    }
    stillUnfound.push_back(pos);
  foundIt:;
    ++fuzzyProgress;
  }
  cerr<<"\r";
  cerr<<(boost::format("Fuzzy found: %|40t|-%10d\n")%fuzzyFound).str();
  fuzzyFound=0;
  unfoundReads.swap(stillUnfound);
  cerr<<(boost::format("Unmatchable reads:%|40t|=%10d (%.2f%%)\n") 
	 % unfoundReads.size() % (100.0*unfoundReads.size()/total)).str();

  cout<<"After: "<<endl;
  /*
  vector<uint32_t> missing{5718050, 1670454,  866724,    1301900,    378768,   2282802,   1728496,  2762187, 985334, 45882};
  for(auto i : missing) {
    cout<<i<<": "<<endl;
    rg.printFastQs(i, fastq);
  }


  exit(1);
  rg.printFastQs(45881, fastq);  
  rg.printFastQs(1728495, fastq);
  */
  /*
  rg.printFastQs(5718000, fastq);  
  FILE *fp=fopen("unfound.fastq", "w");
  BOOST_FOREACH(uint64_t pos, unfoundReads) {
    fastq.seek(pos);
    fastq.getRead(&fqfrag);
    fprintf(fp, "%s\n%s\n", fqfrag.d_nucleotides.c_str(), fqfrag.d_quality.c_str());
  }
  fclose(fp);
  */
  cout<<"After sliding matching: "<<endl;
  printUnmatched(rg, "fuzzy");

  cerr<<"Found "<<locimap.size()<<" varying loci"<<endl;
  uint64_t seriouslyVariable=0;
  boost::format fmt1("%-10d: %3d*%c ");
  string fmt2("                  ");
  int aCount, cCount, tCount, gCount;
  double fraction;
  for(locimap_t::iterator iter = locimap.begin(); iter != locimap.end(); ++iter) {
    unsigned int varcount=variabilityCount(rg, iter->first, iter->second, &fraction);
    if(varcount < 20) 
      continue;
    char c=rg.snippet(iter->first, iter->first+1)[0];
    aCount = cCount = tCount = gCount = 0;
    switch(c) {
    case 'A':
      aCount += rg.d_mapping[iter->first].coverage;
      break;
    case 'C':
      cCount += rg.d_mapping[iter->first].coverage;
      break;
    case 'T':
      tCount += rg.d_mapping[iter->first].coverage;
      break;
    case 'G':
      gCount += rg.d_mapping[iter->first].coverage;
      break;
    }

    cout<< (fmt1 % iter->first % rg.d_mapping[iter->first].coverage % rg.snippet(iter->first, iter->first+1) ).str();
    sort(iter->second.samples.begin(), iter->second.samples.end());

    seriouslyVariable++;
    for(vector<boost::tuple<char,char,char> >::const_iterator j = iter->second.samples.begin(); 
	j != iter->second.samples.end(); ++j) {
      c=j->get<0>();
      switch(c) {
      case 'A':
	aCount++;
	break;
      case 'C':
	cCount++;
	break;
      case 'T':
	tCount++;
	break;
      case 'G':
	gCount++;
	break;
      }
      
      cout<<c;
    }
    cout<<endl<<fmt2;
    for(vector<boost::tuple<char,char,char> >::const_iterator j = iter->second.samples.begin(); 
	j != iter->second.samples.end(); ++j) {
      cout<<j->get<1>();
    }
    cout<<endl<<fmt2;
    for(vector<boost::tuple<char,char,char> >::const_iterator j = iter->second.samples.begin(); 
	j != iter->second.samples.end(); ++j) {
      cout<< (j->get<2>() ? 'R' : '.');
    }

    int tot=iter->second.samples.size() + rg.d_mapping[iter->first].coverage;
    cout<<endl;
    vector<GeneAnnotation> gas=gar.lookup(iter->first);
    if(!gas.empty()) {
      cout<<fmt2<<"Annotation: ";
      BOOST_FOREACH(const GeneAnnotation& ga, gas) {
	cout<<ga.name<<" ["<<ga.tag<<"], ";
      }
      cout<<endl;
    }
    cout<<fmt2<< "Fraction tail: "<<fraction<<", "<<iter->second.samples.size()<<endl;
    cout<<fmt2<< "A: " << aCount*100/tot <<"%, C: "<<cCount*100/tot<<"%, G: "<<gCount*100/tot<<"%, T: "<<tCount*100/tot<<"%"<<endl;

    rg.printFastQs(iter->first, fastq);

  }
  cerr<<"Found "<<seriouslyVariable<<" seriously variable loci"<<endl;
  exit(EXIT_SUCCESS);

}
