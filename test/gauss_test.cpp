#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "gauss.hpp"
#include "gauss_lagrange_2_2d.hpp"
#include "name_reporter.h"
#include <doctest/doctest.h>

const float EPSILON = std::numeric_limits<float>::epsilon();

TEST_CASE("testing gauss") {
  femib::gauss::rule<float, 2> rule =
      femib::gauss::create_gauss_2_2d<float, 2>();

  std::function<float(femib::types::dvec<float, 2>)> f =
      [](femib::types::dvec<float, 2> x) { return 1; };
  float area = femib::gauss::integrate<float, 2>(rule, f);
  CHECK(std::fabs(area - 1.0 / 2.0) <= EPSILON);

  std::function<float(femib::types::dvec<float, 2>)> g =
      [](femib::types::dvec<float, 2> x) { return x(0) + x(1); };
  float integral = femib::gauss::integrate<float, 2>(rule, g);
  CHECK(std::fabs(integral - 1.0 / 3.0) <= EPSILON);
}
