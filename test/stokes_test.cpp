#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../src/affine/affine.hpp"
#include "../src/femib/stokes_steady.hpp"
#include "../src/finite_element/P0_2d1d.hpp"
#include "../src/finite_element/P1+B_2d2d.hpp"
#include "../src/finite_element/P1_2d1d.hpp"
#include "../src/finite_element_space/finite_element_space.hpp"
#include "../src/gauss/gauss.hpp"
#include "../src/gauss/gauss_lagrange_2_2d.hpp"
#include "../src/mesh/mesh.hpp"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <doctest/doctest.h>
#include <iostream>
#include <vector>

#include "utils.hpp"

// TODO is this test still needed? as femib, it doesn't check anything
TEST_CASE("testing femib stokes") {

  femib::gauss::rule<float, 2> rule =
      femib::gauss::create_gauss_2_2d<float, 2>();
  std::string mesh_dir = MESH_DIR;
  femib::types::mesh<float, 2> mesh = femib::mesh::read<float, 2>(
      mesh_dir + "p3.mat", mesh_dir + "t3.mat", mesh_dir + "e3.mat");

  // V
  femib::finite_element::finite_element<float, 2, 2> f_p1_2d2d =
      femib::finite_element::create_finite_element_P1_B_2d2d<float, 2, 2>();
  femib::finite_element_space::finite_element_space<float, 2, 2> v = {f_p1_2d2d,
                                                                      mesh};
  v.nodes = f_p1_2d2d.build_nodes(mesh);

  // Q
  femib::finite_element::finite_element<float, 2, 1> f_p0_2d1d =
      femib::finite_element::create_finite_element_P0_2d1d<float, 2, 1>();
  femib::finite_element_space::finite_element_space<float, 2, 1> q = {f_p0_2d1d,
                                                                      mesh};
  q.nodes = f_p0_2d1d.build_nodes(mesh);

  // STOKES
  femib::stokes_steady::stokes<float, 2> stokes(
      v, q, rule, /*rho=*/1.0f, /*mu=*/1.0f,
      [](const femib::types::dvec<float, 2> &,
         float) -> femib::types::dvec<float, 2> {
        return femib::types::dvec<float, 2>(1.0f, 1.0f);
      });

  {
    std::vector<int> not_edges =
        femib::util::build_not_edges<float, 2, 2>(stokes.V);
    int rowsV_free = static_cast<int>(not_edges.size());
    int rowsQ = static_cast<int>(stokes.Q.nodes.P.size());
    int expected_n = rowsV_free + rowsQ + 1;

    const Eigen::SparseMatrix<float> &A = stokes.solvable_equations.A;
    const Eigen::Matrix<float, Eigen::Dynamic, 1> &b =
        stokes.solvable_equations.b;

    // Expected size: (rowsV_free + rowsQ + 1) -- one Lagrange-multiplier
    // row/column on top of the velocity-reduced, fully-free-pressure system.
    REQUIRE(A.rows() == expected_n);
    REQUIRE(A.cols() == expected_n);
    REQUIRE(b.rows() == expected_n);

    // Symmetric saddle-point structure.
    CHECK((A - Eigen::SparseMatrix<float>(A.transpose())).norm() < 1e-5f);

    // Bottom-right corner is the Tikhonov regularization term, not exactly 0
    CHECK(std::abs(A.coeff(expected_n - 1, expected_n - 1)) < 1e-5f);

    // RHS of "integral(p) = 0" is exactly 0.
    CHECK(b(expected_n - 1) == 0.0f);

    // The augmented row/column (excluding the corner) is non-trivial
    CHECK(Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>(A)
              .block(0, expected_n - 1, expected_n - 1, 1)
              .norm() > 0.0f);
  }

  Eigen::Matrix<float, Eigen::Dynamic, 1> xx =
      femib::stokes_steady::solve<float, 2, 1>(stokes);
  CHECK(xx.allFinite());
  CHECK(xx.topRows(stokes.V.size()).norm() > 0.0f);
}

TEST_CASE("stokes::solve matches curl(sin^2(pi x)sin^2(pi y)) manufactured "
          "solution") {
  // TODO stokes implementation supports only 0-boundary-condition... we need
  //  others?

  const float PI = 3.14159265358979f; // TODO math

  auto u1_exact = [PI](float x, float y) {
    return 2.0f * PI * std::sin(PI * x) * std::sin(PI * x) * std::sin(PI * y) *
           std::cos(PI * y);
  };
  auto u2_exact = [PI](float x, float y) {
    return -2.0f * PI * std::sin(PI * x) * std::sin(PI * y) * std::sin(PI * y) *
           std::cos(PI * x);
  };
  auto p_exact = [PI](float x, float y) {
    return std::sin(PI * x) * std::cos(PI * y);
  };
  auto f1 = [PI](float x, float y) {
    return 0.5 * PI *
           (32.0f * PI * PI * std::sin(PI * x) * std::sin(PI * x) *
                std::sin(PI * y) -
            8.0f * PI * PI * std::sin(PI * y) - std::cos(PI * x)) *
           std::cos(PI * y);
  };
  auto f2 = [PI](float x, float y) {
    return 0.5 * PI *
           (-32.0f * PI * PI * std::sin(PI * y) * std::sin(PI * y) *
                std::cos(PI * x) +
            std::sin(PI * y) + 8.0f * PI * PI * std::cos(PI * x)) *
           std::sin(PI * x);
  };

  int n = 20;
  femib::types::mesh<float, 2> mesh = make_unit_square_mesh(n);

  femib::gauss::rule<float, 2> rule =
      femib::gauss::create_gauss_2_2d<float, 2>();

  femib::finite_element::finite_element<float, 2, 2> f_p1_2d2d =
      femib::finite_element::create_finite_element_P1_B_2d2d<float, 2, 2>();
  femib::finite_element_space::finite_element_space<float, 2, 2> v = {f_p1_2d2d,
                                                                      mesh};
  v.nodes = f_p1_2d2d.build_nodes(mesh);

  femib::finite_element::finite_element<float, 2, 1> f_p1_2d1d =
      femib::finite_element::create_finite_element_P1_2d1d<float, 2, 1>();
  femib::finite_element_space::finite_element_space<float, 2, 1> q = {f_p1_2d1d,
                                                                      mesh};
  q.nodes = f_p1_2d1d.build_nodes(mesh);

  femib::stokes_steady::stokes<float, 2> stokes(
      v, q, rule, /*rho=*/1.0f, /*mu=*/1.0f,
      [f1, f2](const femib::types::dvec<float, 2> &x,
               float) -> femib::types::dvec<float, 2> {
        return femib::types::dvec<float, 2>(f1(x(0), x(1)), f2(x(0), x(1)));
      });

  Eigen::Matrix<float, Eigen::Dynamic, 1> xx =
      femib::stokes_steady::solve<float, 2, 1>(stokes);

  int size_P = mesh.P.size();
  int size_T = mesh.T.size();

  float max_velocity_error = 0.0f;
  for (int i = 0; i < size_P; ++i) {
    float x = mesh.P[i](0);
    float y = mesh.P[i](1);
    float u1_fem = xx(i, 0);
    float u2_fem = xx(size_P + i, 0);
    float u1_ex = u1_exact(x, y);
    float u2_ex = u2_exact(x, y);
    float error_u1 = std::abs(u1_fem - u1_ex);
    float error_u2 = std::abs(u2_fem - u2_ex);
    max_velocity_error = std::max(max_velocity_error, error_u1);
    max_velocity_error = std::max(max_velocity_error, error_u2);
  }
  for (int n_ = 0; n_ < size_T; ++n_) {
    femib::types::dvec<float, 2> c =
        femib::finite_element::find_center_of<float, 2>(mesh[n_]);
    float u1_fem = xx(2 * size_P + n_, 0);
    float u2_fem = xx(2 * size_P + size_T + n_, 0);
    float u1_ex = u1_exact(c(0), c(1));
    float u2_ex = u2_exact(c(0), c(1));
    float error_u1 = std::abs(u1_fem - u1_ex);
    float error_u2 = std::abs(u2_fem - u2_ex);
    max_velocity_error = std::max(max_velocity_error, error_u1);
    max_velocity_error = std::max(max_velocity_error, error_u2);
  }
  CHECK(max_velocity_error < 0.15f);

  // Stokes pressure is only defined up to an additive constant for
  // pure-Dirichlet-velocity boundary conditions, so subtract the mean from
  // both the computed and the exact pressure before comparing.
  int rows_v = stokes.V.nodes.P.size(); // combined velocity DOF count
  int n_p = q.nodes.P.size();
  float mean_computed_p = 0.0f;
  float mean_exact_p = 0.0f;
  for (int n_ = 0; n_ < n_p; ++n_) {
    mean_computed_p += xx(rows_v + n_, 0);
    mean_exact_p += p_exact(q.nodes.P[n_](0), q.nodes.P[n_](1));
  }
  mean_computed_p /= n_p;
  mean_exact_p /= n_p;

  float max_pressure_error = 0.0f;
  for (int n_ = 0; n_ < n_p; ++n_) {
    float raw = xx(rows_v + n_, 0);
    float computed = raw - mean_computed_p;
    float exact = p_exact(q.nodes.P[n_](0), q.nodes.P[n_](1)) - mean_exact_p;
    float err = std::abs(computed - exact);
    max_pressure_error = std::max(max_pressure_error, err);
  }

  CHECK(max_pressure_error < 6.0f);
}
