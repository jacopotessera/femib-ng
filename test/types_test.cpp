#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "name_reporter.h"
#include "types.hpp"
#include <doctest/doctest.h>

TEST_CASE("testing types") {
  femib::types::dvec<float, 2> v = {1.0, 2.0};
  CHECK(v(0) == doctest::Approx(1.0));
  CHECK(v(1) == doctest::Approx(2.0));

  femib::types::F<float, 2, 1> f;
  f.x = [](const femib::types::dvec<float, 2> &x) {
    return femib::types::dvec<float, 1>{x(0) + x(1)};
  };
  femib::types::dvec<float, 2> p = {2.0, 3.0};
  CHECK(f.x(p)(0) == doctest::Approx(5.0));
}

TEST_CASE("testing mesh init") {
  femib::types::mesh<float, 2> mesh;
  mesh.P = {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}};
  mesh.T = {{0, 1, 2}};
  mesh.E = {};

  mesh.init();
  size_t n_after_first = mesh.N.size();
  mesh.init();
  size_t n_after_second = mesh.N.size();

  CHECK(n_after_first > 0);
  CHECK(n_after_second == n_after_first);
}
