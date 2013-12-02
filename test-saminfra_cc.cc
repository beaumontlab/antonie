#include <boost/test/unit_test.hpp>
#include "saminfra.hh"
BOOST_AUTO_TEST_SUITE(saminfra_hh)
using std::string;

BOOST_AUTO_TEST_CASE(test_bamCompress) {
  BOOST_CHECK_EQUAL(bamCompress("AAAA"), string("\x11\x11", 2));
  BOOST_CHECK_EQUAL(bamCompress("CCCC"), string("\x22\x22", 2));
  BOOST_CHECK_EQUAL(bamCompress("ACACACAC"), string("\x12\x12\x12\x12", 4));
  BOOST_CHECK_EQUAL(bamCompress("NNNN"), string("\xff\xff", 2));
  BOOST_CHECK_EQUAL(bamCompress("PPPP"), string("\xff\xff", 2));
}


BOOST_AUTO_TEST_SUITE_END()
