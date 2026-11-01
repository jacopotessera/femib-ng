#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../src/affine/affine.hpp"
#include "../src/femib/femib.hpp"
#include "../src/femib/poisson.hpp"
#include "../src/finite_element/P1_2d1d.hpp"
#include "../src/finite_element_space/finite_element_space.hpp"
#include "../src/gauss/gauss.hpp"
#include "../src/gauss/gauss_2_2d.hpp"
#include "../src/mesh/mesh.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <cmath>
#include <doctest/doctest.h>
#include <numbers>
#include <vector>

#include "utils.hpp"

TEST_CASE("testing femib poisson") {

  femib::gauss::rule<float, 2> rule =
      femib::gauss::create_gauss_2_2d<float, 2>();
  std::string mesh_dir = MESH_DIR;
  femib::types::mesh<float, 2> mesh = femib::mesh::read<float, 2>(
      mesh_dir + "p5.mat", mesh_dir + "t5.mat", mesh_dir + "e5.mat");

  femib::finite_element::finite_element<float, 2, 1> f =
      femib::finite_element::create_finite_element_P1_2d1d<float, 2, 1>();

  femib::finite_element_space::finite_element_space<float, 2, 1> s(f, mesh);

  femib::poisson::poisson<float, 2, 1> poisson(s, rule);

  Eigen::Matrix<float, Eigen::Dynamic, 1> xx =
      femib::poisson::solve<float, 2, 1>(poisson);

  CHECK_UNARY(xx.allFinite());
  // TODO check? it is a known problem with known solutions. there is a test
  //  below. is this test still needed?
}

TEST_CASE("testing find_not_edge_nodes") {
  // 4-triangle square mesh: 4 boundary corner nodes (0,1,2,3) plus one
  // interior center node (4).
  femib::types::mesh<float, 2> mesh(
      {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}, {0.5, 0.5}},
      {{0, 1, 4}, {1, 2, 4}, {2, 3, 4}, {3, 0, 4}}, {0, 1, 2, 3});

  femib::finite_element::finite_element<float, 2, 1> f =
      femib::finite_element::create_finite_element_P1_2d1d<float, 2, 1>();

  femib::finite_element_space::finite_element_space<float, 2, 1> s(f, mesh);

  std::vector<int> not_edge_nodes =
      femib::util::find_not_edge_nodes<float, 2, 1>(s);
  CHECK_EQ(not_edge_nodes.size() + s.nodes.E.size(), s.nodes.P.size());
  for (int i : not_edge_nodes) {
    CHECK_EQ(std::find(s.nodes.E.begin(), s.nodes.E.end(), i), s.nodes.E.end());
  }
}

TEST_CASE("testing triplets2dense") {
  std::vector<Eigen::Triplet<float>> triplets = {
      {0, 0, 2.0f}, {1, 2, 3.5f}, {2, 1, -1.0f}};

  Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> dense =
      femib::util::triplets2dense<float>(triplets, 3, 3);

  REQUIRE_EQ(dense.rows(), 3);
  REQUIRE_EQ(dense.cols(), 3);
  CHECK_EQ(dense(0, 0), doctest::Approx(2.0f));
  CHECK_EQ(dense(1, 2), doctest::Approx(3.5f));
  CHECK_EQ(dense(2, 1), doctest::Approx(-1.0f));

  for (int i = 0; i < dense.rows(); ++i) {
    for (int j = 0; j < dense.cols(); ++j) {
      bool listed = (i == 0 && j == 0) || (i == 1 && j == 2) ||
                    (i == 2 && j == 1); // TODO eh. use find_if?
      if (!listed) {
        CHECK_EQ(dense(i, j), doctest::Approx(0.0f));
      }
    }
  }
}

TEST_CASE("tesint triplets2sparse") {
  std::vector<Eigen::Triplet<float>> triplets = {
      {0, 0, 2.0f}, {1, 2, 3.5f}, {2, 1, -1.0f}};

  Eigen::SparseMatrix<float> sparse =
      femib::util::triplets2sparse<float>(triplets, 3, 3);

  REQUIRE_EQ(sparse.rows(), 3);
  REQUIRE_EQ(sparse.cols(), 3);
  REQUIRE_EQ(sparse.nonZeros(), 3);
  CHECK_EQ(sparse.coeff(0, 0), doctest::Approx(2.0f));
  CHECK_EQ(sparse.coeff(1, 2), doctest::Approx(3.5f));
  CHECK_EQ(sparse.coeff(2, 1), doctest::Approx(-1.0f));

  for (int i = 0; i < sparse.rows(); ++i) {
    for (int j = 0; j < sparse.cols(); ++j) {
      bool listed = (i == 0 && j == 0) || (i == 1 && j == 2) ||
                    (i == 2 && j == 1); // TODO as above, use find or similar
      if (!listed) {
        CHECK_EQ(sparse.coeff(i, j), doctest::Approx(0.0f));
      }
    }
  }

  // same result as triplets2dense
  Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> dense_check =
      femib::util::triplets2dense<float>(triplets, 3, 3);
  Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> sparse_densified =
      Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>(sparse);
  CHECK_EQ((dense_check - sparse_densified).norm(), doctest::Approx(0.0f));
}

TEST_CASE("poisson::solve matches sin(pi x)sin(pi y) manufactured solution") {
  // Manufactured solution on the unit square with homogeneous Dirichlet BC
  // (automatically satisfied: sin(pi*0) = sin(pi*1) = 0 on every edge):
  //   u(x,y)      = sin(pi x) sin(pi y)
  //   -Delta(u)   = 2 pi^2 sin(pi x) sin(pi y)
  // (standard, single-term textbook derivation: each second partial
  // derivative contributes -pi^2 sin(pi x) sin(pi y), so
  // -Delta(u) = -(u_xx + u_yy) = 2 pi^2 sin(pi x) sin(pi y).)
  const float PI = std::numbers::pi_v<float>;
  auto u_exact = [PI](const femib::types::dvec<float, 2> &x) {
    return std::sin(PI * x(0)) * std::sin(PI * x(1));
  };
  auto forcing = [PI](const femib::types::dvec<float, 2> &x) {
    return 2.0f * PI * PI * std::sin(PI * x(0)) * std::sin(PI * x(1));
  };

  int n = 16;
  femib::types::mesh<float, 2> mesh = make_unit_square_mesh(n);

  femib::gauss::rule<float, 2> rule =
      femib::gauss::create_gauss_2_2d<float, 2>();

  femib::finite_element::finite_element<float, 2, 1> fe =
      femib::finite_element::create_finite_element_P1_2d1d<float, 2, 1>();

  femib::finite_element_space::finite_element_space<float, 2, 1> s(fe, mesh);

  // the manufactured forcing tested against each basis function (e=1)
  auto ggg = [forcing](femib::types::function_with_gradient<float, 2, 1> a) {
    return [a, forcing](const femib::types::dvec<float, 2> &x) {
      return forcing(x) * a.value(x)(0);
    };
  };

  femib::poisson::poisson<float, 2, 1> problem(
      s, rule, femib::poisson::default_boundary<float, 2>(), ggg);

  Eigen::Matrix<float, Eigen::Dynamic, 1> xx =
      femib::poisson::solve<float, 2, 1>(problem);

  float max_error = 0.0f;
  for (int i = 0; i < std::ssize(s.nodes.P); ++i) {
    float exact = u_exact(s.nodes.P[i]);
    float error = std::abs(xx(i, 0) - exact);
    max_error = std::max(max_error, error);
  }
  CHECK_LT(max_error, 0.05f);
}
