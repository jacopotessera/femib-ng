#define BOOST_TEST_MODULE decorator_10
#include <boost/test/unit_test.hpp>

#include "affine.hpp"

BOOST_AUTO_TEST_SUITE(suite1,
  * boost::unit_test::expected_failures(1))

BOOST_AUTO_TEST_CASE(FailTest,
    * boost::unit_test::expected_failures(1))
{
    dvec p1;
    p1 << 1,2;
    dvec p2;
    p2 << 3,0;
    dvec p3;
    p3 << 2,4;

    dvec x;
    x << 1,0;
    
    std::vector<dvec> t;
    t.push_back(p1);
    t.push_back(p2);
    t.push_back(p3);
    dvec y = affine(t, x);

    BOOST_CHECK_EQUAL(y(0), 12);
}

BOOST_AUTO_TEST_CASE(PassTest)
{
    dvec p1;
    p1 << 1,2;
    dvec p2;
    p2 << 3,0;
    dvec p3;
    p3 << 2,4;

    dvec x;
    x << 1,0;

    std::vector<dvec> t;
    t.push_back(p1);
    t.push_back(p2);
    t.push_back(p3);
    dvec y = affine(t, x);

    BOOST_CHECK_EQUAL(y(0), 3);
}

BOOST_AUTO_TEST_SUITE_END()
