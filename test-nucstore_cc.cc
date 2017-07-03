#include <boost/test/unit_test.hpp>
#include "dnamisc.hh"
#include "nucstore.hh"

BOOST_AUTO_TEST_SUITE(nucstore_hh)

BOOST_AUTO_TEST_CASE(test_nucstore) {
  NucleotideStore ns;
  ns.append('A');
  BOOST_CHECK_EQUAL(ns.get(0),'A');
  ns.append('C');
  BOOST_CHECK_EQUAL(ns.get(0),'A');
  BOOST_CHECK_EQUAL(ns.get(1),'C');
  ns.append('G');
  BOOST_CHECK_EQUAL(ns.get(2),'G');
  ns.append('T');

  BOOST_CHECK_EQUAL(ns.get(3),'T');

  BOOST_CHECK_EQUAL(ns.get(0),'A');
  ns.append('C');
  BOOST_CHECK_EQUAL(ns.get(0),'A');
  BOOST_CHECK_EQUAL(ns.get(1),'C');
  ns.append('G');
  BOOST_CHECK_EQUAL(ns.get(2),'G');
  ns.append('T');

  BOOST_CHECK_EQUAL(ns.get(3),'T');


  BOOST_CHECK_EQUAL(ns.size(), 7);
  
  // AACG TAACG

  NucleotideStore sep;
  sep.append("ACGTCGT");
  BOOST_CHECK_EQUAL(ns, sep);

  
  
}

BOOST_AUTO_TEST_SUITE_END()
