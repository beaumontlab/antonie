#include <boost/test/unit_test.hpp>
#include "dnamisc.hh"
BOOST_AUTO_TEST_SUITE(misc_hh)

BOOST_AUTO_TEST_CASE(test_kmerMapper) {
  BOOST_CHECK_EQUAL(kmerMapper("AAAA", 0, 4), 0U);
  BOOST_CHECK_EQUAL(kmerMapper("AAAAAAAA", 0, 8), 0U);
  BOOST_CHECK_EQUAL(kmerMapper("AAAAAAAAAAAA", 0, 12), 0U);
  BOOST_CHECK_EQUAL(kmerMapper("AAAAAAAAAAAAAAAA", 0, 16), 0U);

  BOOST_CHECK_EQUAL(kmerMapper("CCCC", 0, 4), 85U);
  BOOST_CHECK_EQUAL(kmerMapper("CCCCCCCC", 0, 8), 21845U);
  BOOST_CHECK_EQUAL(kmerMapper("CCCCCCCCCCCC", 0, 12), 5592405U);
  BOOST_CHECK_EQUAL(kmerMapper("CCCCCCCCCCCCCCCC", 0, 16), 1431655765U);

  BOOST_CHECK_EQUAL(DNAToAminoAcid("GCC"), 'A');
}

BOOST_AUTO_TEST_SUITE_END()
