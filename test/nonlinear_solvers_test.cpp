#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../src/femib/nonlinear_solvers.hpp"
#include <Eigen/Dense>
#include <doctest/doctest.h>
#include <stdexcept>

using vector_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;

TEST_CASE("picard_solver converges to the fixed point of cos") {
  femib::util::picard_solver<double> solver(100, 1e-12);
  vector_t start = vector_t::Constant(3, 1.0);
  vector_t x = solver.solve(
      start,
      [](const vector_t &u) -> vector_t { return u.array().cos().matrix(); },
      [](const vector_t &full) { return full; });
  for (int i = 0; i < 3; ++i) {
    CHECK_EQ(x(i), doctest::Approx(0.7390851332151607).epsilon(1e-9));
  }
}

TEST_CASE("picard_solver never stops on the first step") {
  femib::util::picard_solver<double> solver(10, 1e-6);
  int calls = 0;
  vector_t start = vector_t::Constant(2, 5.0);
  vector_t x = solver.solve(
      start,
      [&calls](const vector_t &u) -> vector_t {
        ++calls;
        return u;
      },
      [](const vector_t &full) { return full; });
  CHECK_EQ(calls, 2);
  CHECK_EQ((x - start).norm(), 0.0);
}

TEST_CASE("picard_solver stops after max_iters steps") {
  femib::util::picard_solver<double> solver(3, 1e-6);
  int calls = 0;
  vector_t start = vector_t::Zero(2);
  vector_t x = solver.solve(
      start,
      [&calls](const vector_t &u) -> vector_t {
        ++calls;
        return (u.array() + 1.0).matrix();
      },
      [](const vector_t &full) { return full; });
  CHECK_EQ(calls, 3);
  CHECK_EQ((x - vector_t::Constant(2, 3.0)).norm(), 0.0);
}

TEST_CASE("picard_solver rejects max_iters <= 0") {
  femib::util::picard_solver<double> solver(0, 1e-6);
  vector_t start = vector_t::Zero(1);
  auto identity = [](const vector_t &u) { return u; };
  CHECK_THROWS_AS(solver.solve(start, identity, identity),
                  std::invalid_argument);
}
