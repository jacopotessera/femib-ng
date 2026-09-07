#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../src/affine/affine.hpp"
#include "../src/femib/femib.hpp"
#include "../src/femib/poisson.hpp"
#include "../src/finite_element/P1_2d1d.hpp"
#include "../src/finite_element_space/finite_element_space.hpp"
#include "../src/gauss/gauss.hpp"
#include "../src/gauss/gauss_lagrange_2_2d.hpp"
#include "../src/mesh/mesh.hpp"
#include "../src/read/read.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <cmath>
#include <doctest/doctest.h>
#include <iostream>
#include <vector>

#include "utils.hpp"

TEST_CASE("testing femib poisson") {

  femib::gauss::rule<float, 2> rule =
      femib::gauss::create_gauss_2_2d<float, 2>();
  std::string mesh_dir = MESH_DIR;
  femib::types::mesh<float, 2> mesh = femib::mesh::read<float, 2>(
      mesh_dir + "p5.mat", mesh_dir + "t5.mat", mesh_dir + "e5.mat");
  mesh.init();

  femib::finite_element::finite_element<float, 2, 1> f =
      femib::finite_element::create_finite_element_P1_2d1d<float, 2, 1>();

  femib::finite_element_space::finite_element_space<float, 2, 1> s = {
      .finite_element = f, .mesh = mesh};
  s.nodes = f.build_nodes(mesh); // TODO do this in constructor?

  femib::poisson::poisson<float, 2, 1> poisson = {s};

  femib::poisson::init<float, 2, 1>(poisson, rule);
  auto end = std::chrono::steady_clock::now();

  Eigen::Matrix<float, Eigen::Dynamic, 1> xx =
      femib::poisson::solve<float, 2, 1>(poisson);

  CHECK(xx.allFinite());
  // TODO check? it is a known problem with known solutions. there is a test
  //  below. is this test still needed?
}

TEST_CASE("testing build_not_edges") {
  // 4-triangle square mesh: 4 boundary corner nodes (0,1,2,3) plus one
  // interior center node (4).
  femib::types::mesh<float, 2> mesh = {
      .P = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}, {0.5, 0.5}},
      .T = {{0, 1, 4}, {1, 2, 4}, {2, 3, 4}, {3, 0, 4}},
      .E = {0, 1, 2, 3}};
  mesh.init();

  femib::finite_element::finite_element<float, 2, 1> f =
      femib::finite_element::create_finite_element_P1_2d1d<float, 2, 1>();

  femib::finite_element_space::finite_element_space<float, 2, 1> s = {f, mesh};
  s.nodes = f.build_nodes(mesh);

  std::vector<int> not_edges = femib::util::build_not_edges<float, 2, 1>(s);
  // TODO https://github.com/doctest/doctest/issues/1004
  CHECK(not_edges.size() + s.nodes.E.size() == s.nodes.P.size());
  for (int i : not_edges) {
    CHECK(std::find(s.nodes.E.begin(), s.nodes.E.end(), i) == s.nodes.E.end());
  }
}

TEST_CASE("testing triplets2dense") {
  std::vector<Eigen::Triplet<float>> triplets = {
      {0, 0, 2.0f}, {1, 2, 3.5f}, {2, 1, -1.0f}};

  Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> dense =
      femib::util::triplets2dense<float>(triplets, 3, 3);

  REQUIRE(dense.rows() == 3);
  REQUIRE(dense.cols() == 3);
  CHECK(dense(0, 0) == doctest::Approx(2.0f));
  CHECK(dense(1, 2) == doctest::Approx(3.5f));
  CHECK(dense(2, 1) == doctest::Approx(-1.0f));

  for (int i = 0; i < dense.rows(); ++i) {
    for (int j = 0; j < dense.cols(); ++j) {
      bool listed = (i == 0 && j == 0) || (i == 1 && j == 2) ||
                    (i == 2 && j == 1); // TODO eh. use find_if?
      if (!listed) {
        CHECK(dense(i, j) == doctest::Approx(0.0f));
      }
    }
  }
}

TEST_CASE("poisson::solve matches sin(pi x)sin(pi y) manufactured solution") {
  // Manufactured solution on the unit square with homogeneous Dirichlet BC
  // (automatically satisfied: sin(pi*0) = sin(pi*1) = 0 on every edge):
  //   u(x,y)      = sin(pi x) sin(pi y)
  //   -Delta(u)   = 2 pi^2 sin(pi x) sin(pi y)
  // (standard, single-term textbook derivation: each second partial
  // derivative contributes -pi^2 sin(pi x) sin(pi y), so
  // -Delta(u) = -(u_xx + u_yy) = 2 pi^2 sin(pi x) sin(pi y).)
  const float PI = 3.14159265358979f; // TODO math?
  auto u_exact = [PI](const femib::types::dvec<float, 2> &x) {
    return std::sin(PI * x(0)) * std::sin(PI * x(1));
  };
  auto forcing = [PI](const femib::types::dvec<float, 2> &x) {
    return 2.0f * PI * PI * std::sin(PI * x(0)) * std::sin(PI * x(1));
  };

  int n = 16;
  femib::types::mesh<float, 2> mesh = make_unit_square_mesh(n);
  mesh.init();

  femib::gauss::rule<float, 2> rule =
      femib::gauss::create_gauss_2_2d<float, 2>();

  femib::finite_element::finite_element<float, 2, 1> fe =
      femib::finite_element::create_finite_element_P1_2d1d<float, 2, 1>();

  femib::finite_element_space::finite_element_space<float, 2, 1> s = {fe, mesh};
  s.nodes = fe.build_nodes(mesh);

  // TODO this note says: poisson interface is wrong...
  // NOTE on how this test wires in the forcing term: femib::poisson::init's
  // existing pipeline calls build_diagonal with femib::poisson::external_force
  // itself as the RHS functor, and build_diagonal always invokes that functor
  // on a basis function of the SAME finite element space used for the
  // stiffness matrix -- there is no hook to inject an independent, spatially
  // varying source field through that path (it always assembles integral(phi_i)
  // dx, i.e. an implicit constant unit forcing f===1). To genuinely exercise
  // this manufactured solution's nonzero, non-constant forcing term -- and
  // specifically the out-of-bounds `external_force` code path this task
  // fixes -- this test assembles the RHS itself: it still calls
  // femib::poisson::external_force explicitly (so the fix under test is on
  // the exact path exercised) and multiplies its result by the manufactured
  // `forcing` field, then hands the assembled dM/dF/dB to femib::poisson::solve
  // exactly as femib::poisson::init would have.
  auto ggg = [forcing](femib::types::F<float, 2, 1> a) {
    return [a, forcing](const femib::types::dvec<float, 2> &x) {
      return forcing(x) * femib::poisson::external_force<float, 2, 1>(a)(x);
    };
  };

  femib::util::build_diagonal_result<float> result =
      femib::util::build_diagonal<float, 2, 1>(
          s, rule, femib::poisson::ddot<float, 2, 1>, ggg);

  std::function<float(femib::types::dvec<float, 2>)> zero_boundary =
      [](const femib::types::dvec<float, 2> &) { return 0.0f; };
  std::vector<Eigen::Triplet<float>> B =
      femib::util::build_edges<float, 2, 1>(s, zero_boundary);
  std::vector<int> not_edges = femib::util::build_not_edges<float, 2, 1>(s);

  femib::poisson::poisson<float, 2, 1> problem;
  problem.V = s;
  problem.dB = femib::util::triplets2dense<float>(B, s.nodes.P.size(), 1);
  problem.dM = femib::util::triplets2dense<float>(result.M, s.nodes.P.size(),
                                                  s.nodes.P.size());
  problem.dF =
      femib::util::triplets2dense<float>(result.F, s.nodes.P.size(), 1);

  Eigen::Matrix<float, Eigen::Dynamic, 1> xx =
      femib::poisson::solve<float, 2, 1>(problem);

  float max_error = 0.0f;
  for (int i = 0; i < s.nodes.P.size(); ++i) {
    float exact = u_exact(s.nodes.P[i]);
    float error = std::abs(xx(i, 0) - exact);
    max_error = std::max(max_error, error);
  }
  CHECK(max_error < 0.05f);
}
