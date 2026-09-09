#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../src/types/differential_operation.hpp"
#include <doctest/doctest.h>

TEST_CASE("convection((w), (u)) matches (w . grad) u for a known polynomial "
          "field, evaluated by hand") {
  femib::types::F<float, 2, 2> w;
  w.x = [](const femib::types::dvec<float, 2> &) {
    return femib::types::dvec<float, 2>(2.0f, 3.0f);
  };
  w.dx = [](const femib::types::dvec<float, 2> &)
      -> femib::types::rmat<float, 2, 2> {
    return femib::types::rmat<float, 2, 2>::Zero();
  };

  femib::types::F<float, 2, 2> u;
  u.x = [](const femib::types::dvec<float, 2> &x) {
    return femib::types::dvec<float, 2>(x(0) * x(0), x(0) * x(1));
  };
  u.dx = [](const femib::types::dvec<float, 2> &x)
      -> femib::types::rmat<float, 2, 2> {
    femib::types::rmat<float, 2, 2> J;
    J << 2 * x(0), x(1), 0, x(0);
    return J;
  };

  auto conv = convection<float, 2>(w, u);
  femib::types::dvec<float, 2> at(1.5f, -0.5f);
  femib::types::dvec<float, 2> result = conv(at);

  CHECK_LT(std::abs(result(0) - 4.0f * at(0)), 1e-4f);
  CHECK_LT(std::abs(result(1) - (2.0f * at(1) + 3.0f * at(0))), 1e-4f);
}
