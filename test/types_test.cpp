#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "name_reporter.h"
#include "types.hpp"
#include <doctest/doctest.h>

TEST_CASE("testing types") {
  femib::types::dvec<float, 2> v = {1.0, 2.0};
  CHECK_EQ(v(0), doctest::Approx(1.0));
  CHECK_EQ(v(1), doctest::Approx(2.0));

  femib::types::function_with_gradient<float, 2, 1> f;
  f.value = [](const femib::types::dvec<float, 2> &x) {
    return femib::types::dvec<float, 1>{x(0) + x(1)};
  };
  femib::types::dvec<float, 2> p = {2.0, 3.0};
  CHECK_EQ(f.value(p)(0), doctest::Approx(5.0));
}

TEST_CASE("testing mesh constructor populates N") {
  femib::types::mesh<float, 2> mesh({{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}},
                                    {{0, 1, 2}}, {});

  CHECK_EQ(mesh.N.size(), mesh.T.size());
}
