#include <boost/test/unit_test.hpp>
#include "misc.hh"
BOOST_AUTO_TEST_SUITE(misc_hh)

BOOST_AUTO_TEST_CASE(test_VarMeanEstimator) {
	VarMeanEstimator vme;
	vme(0);
	BOOST_CHECK_CLOSE(mean(vme), 0.0, 0.001);
	BOOST_CHECK_CLOSE(variance(vme), 0.0, 0.001);

	for(auto d : {1,2,3,4})
		vme(d);
	BOOST_CHECK_CLOSE(mean(vme), 2.0, 0.001);
	BOOST_CHECK_CLOSE(variance(vme), 2.0, 0.001);

}
BOOST_AUTO_TEST_CASE(test_reverseNucleotides) {
	std::string tst{"TTTTGGGCCA"};
	reverseNucleotides(&tst);
	BOOST_CHECK_EQUAL(tst, "TGGCCCAAAA");
	tst.clear();
	reverseNucleotides(&tst);
	BOOST_CHECK_EQUAL(tst, "");
	
}

BOOST_AUTO_TEST_SUITE_END()
