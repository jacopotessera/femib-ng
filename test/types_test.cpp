#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "name_reporter.h"
#include "operators.hpp"
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

TEST_CASE("gradient_double_dot sums every entry of the two gradients") {
  using vec = femib::types::dvec<double, 2>;
  auto constant_gradient = [](femib::types::rmat<double, 2, 2> g) {
    femib::types::function_with_gradient<double, 2, 2> f;
    f.value = [](const vec &) { return vec::Zero().eval(); };
    f.gradient = [g](const vec &) { return g; };
    return f;
  };
  femib::types::rmat<double, 2, 2> ga, gb;
  ga << 1, 2, 3, 4;
  gb << 5, 6, 7, 8;
  auto vector_product = femib::operators::gradient_double_dot<double, 2, 2>(
      constant_gradient(ga), constant_gradient(gb));
  CHECK_EQ(vector_product(vec(0.3, 0.7)), doctest::Approx(70.0));

  femib::types::function_with_gradient<double, 2, 1> a, b;
  a.value = b.value = [](const vec &) {
    return femib::types::dvec<double, 1>::Zero().eval();
  };
  a.gradient = [](const vec &) { return vec(1, 2); };
  b.gradient = [](const vec &) { return vec(3, 4); };
  auto scalar_product =
      femib::operators::gradient_double_dot<double, 2, 1>(a, b);
  CHECK_EQ(scalar_product(vec(0.3, 0.7)), doctest::Approx(11.0));
}
