#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../src/gauss/gauss.hpp"
#include <doctest/doctest.h>

TEST_CASE("testing gauss") {
  femib::gauss::node<float, 2> node1 = {1.0 / 6.0, {1.0 / 6.0, 1.0 / 6.0}};
  femib::gauss::node<float, 2> node2 = {1.0 / 6.0, {1.0 / 6.0, 2.0 / 3.0}};
  femib::gauss::node<float, 2> node3 = {1.0 / 6.0, {2.0 / 3.0, 1.0 / 6.0}};
  femib::gauss::rule<float, 2> rule = {{node1, node2, node3}};

  std::function<float(femib::types::dvec<float, 2>)> f =
      [](femib::types::dvec<float, 2> x) { return 1; };
  float area = femib::gauss::integrate<float, 2>(rule, f);

  CHECK(area == 1.0 / 2.0);
}
