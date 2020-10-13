#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include "affine.hpp"

TEST_CASE("testing affine") {
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

    CHECK(y(0) == 3);
    CHECK(y(0) == 4);
}

