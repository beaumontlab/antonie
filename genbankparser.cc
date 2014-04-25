#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include "geneannotated.hh"
using std::cout;
using std::endl;

namespace {

  struct State
  {
    void clear()
    {
      startLocus = stopLocus = 0;
      strand = true;
      features.clear();
    }
    uint32_t startLocus;
    uint32_t stopLocus;
    bool strand;
    std::string kind;
    std::vector<std::pair<std::string, std::string> > features;
  } state;

  std::vector<GeneAnnotation> g_ret;

  void reportKind(const std::string& kind)
  {
    if(!state.kind.empty()) {
      GeneAnnotation ga;
      ga.startPos = state.startLocus;
      ga.stopPos = state.stopLocus;
      ga.strand = state.strand;
      ga.type=state.kind;
      ga.gene=false;
      /*
	cout<<"Should emit '"<<state.kind<<"', "<<state.startLocus<<" - " << state.stopLocus<< " "<< (state.strand ? '+' : '-')<<endl;
      */
      for(auto& a: state.features) {
	if(a.first!="translation") {
	  boost::replace_all(a.second, "\n", " ");
	  boost::replace_all(a.second, "'", "\\'");
	  ga.tag+=a.second;
	  ga.tag+=+", ";
	}
      }
      if(ga.type=="CDS" || ga.type=="gene")
	ga.gene=true;

      if(ga.type!="source")
	g_ret.push_back(ga);
      state.clear();
    }
    state.kind=kind;
  }

  // XXX FIXME, we get a lot of these for order() and join()!
  void startLocus(int start)
  {
    state.startLocus = start;

  }

  void stopLocus(int stop)
  {
    state.stopLocus = stop;

  }

  void complement()
  {
    state.strand=false;
  }


  void variable(const std::string& var)
  {
    state.features.push_back({var, std::string()});
  }

  void value(int val)
  {
    state.features.rbegin()->second=std::to_string(val);
  }

  void stringValue(const std::string& val)
  {
    state.features.rbegin()->second=val;
  }


}

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;


std::vector<GeneAnnotation> parseGenBankString(const std::string& bank)
{
    auto first = bank.begin();
    auto last=bank.end();
    g_ret.clear();

    using qi::phrase_parse;
    using qi::lit;
    using qi::lexeme;
    using qi::alpha;
    using qi::char_;
    using qi::int_;
    using ascii::space;

    qi::rule<std::string::const_iterator, std::string(), ascii::space_type> quoted_string, unquoted_string, number_range;
    quoted_string %= lexeme['"' >> +(char_ - '"') >> '"'];
    unquoted_string %= lexeme[+(alpha | char_('_'))];
    number_range %= lexeme[-char_('<') >> int_[startLocus] >> lit("..") >> -char_('>') >> int_[stopLocus]];


    auto base_range = (number_range) |
      (lit("complement(")[complement] >> number_range >> char_(')')) |
      (lit("order(") >> *(number_range >> -char_(',') ) >> lit(")")  ) |
      (lit("join(") >> *(number_range >> -char_(',') ) >> lit(")")  ) |
      (lit("complement(order(")[complement] >> *(number_range >> -char_(',') ) >> lit("))")  ) |
      (lit("complement(join(")[complement] >> *(number_range >> -char_(',') ) >> lit("))"));


    // anti_codon needs to be treated like transl_except 
    bool r = phrase_parse(
	first,                         
	last,                          
	*((unquoted_string[reportKind] >> 
	   base_range
	   >> *(char_('/') >> (
			       (lit("transl_except=(pos:") >> base_range >> char_(',') >> lit("aa:") >> unquoted_string >> lit(")"))  |
			       (unquoted_string[variable] >> -(char_('=') >> (int_[value] | (quoted_string[stringValue] )  ))) 

			       )
		)
	   )
	  ) 
	
	
	,
	space                           /*< the skip-parser >*/
			  );
    if (!r || first != last) {// fail if we did not get a full match	
      //      cout<<"Failed at: '"<<bank.substr(last-first, 100)<<"'"<<endl;
      throw std::runtime_error("Failed to parse genbank string at byte "+std::to_string(last-first) +" of " +std::to_string(bank.size()));
    }
    reportKind("");

    return g_ret;
}

