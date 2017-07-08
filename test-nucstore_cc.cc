#include <boost/test/unit_test.hpp>
#include "dnamisc.hh"
#include "nucstore.hh"
#include <iostream>

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


  sep.set(0, 'C');
  BOOST_CHECK_EQUAL(sep.get(0), 'C');


  sep.set(4, 'C');
  BOOST_CHECK_EQUAL(sep.get(0), 'C');

  
  
}


BOOST_AUTO_TEST_CASE(test_delta) {
  using namespace std;
  NucleotideStore a("ACGTTGCA"), b("ACGTTTCA"), c;
  auto ds = a.getDelta(b);
  cout<<a<<endl<<b<<endl;
  for(const auto& d : ds) {
    cout<<d<<endl;
  }
  cout<<"---"<<endl;

  vector<NucleotideStore::Delta> expected({{(uint32_t)5, 'T', NucleotideStore::Delta::Action::Replace}});
  BOOST_CHECK(ds==expected);
  
  auto ds2= a.getDelta(c);
  for(const auto& d : ds2) {
    cout<<d<<endl;
  }
  
  NucleotideStore d("AGCCTTTCCGGA"), e("AGCCTTTCCCGGA");
  auto ds3 = d.getDelta(e);
  cout<<"---"<<endl;

  for(const auto& d : ds3) {
    cout<<d<<endl;
  }
  
  NucleotideStore f("AGCCTTTCCCGGGA");
  cout<<"---"<<endl;
  for(const auto& de : d.getDelta(f))
    cout<<de<<endl;

  
}
BOOST_AUTO_TEST_SUITE_END()
