#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../src/finite_element/P1+B_2d2d.hpp"
#include "../src/finite_element_space/finite_element_space.hpp"
#include "../src/ib/coupling.hpp"
#include "../src/mesh/mesh.hpp"
#include <doctest/doctest.h>

#include "structure.hpp"
#include "utils.hpp"

TEST_CASE("interpolate_velocity exactly reproduces a linear velocity field") {
  // MINI element reproduce any affine field exactly
  femib::types::mesh<float, 2> mesh = make_unit_square_mesh(6);
  mesh.init();

  femib::finite_element::finite_element<float, 2, 2> fe =
      femib::finite_element::create_finite_element_P1_B_2d2d<float, 2, 2>();
  femib::finite_element_space::finite_element_space V = {.finite_element = fe,
                                                         .mesh = mesh};
  V.nodes = fe.build_nodes(mesh);

  auto u_exact = [](const float x, const float y) {
    return femib::types::dvec<float, 2>(1.0f + 2.0f * x - 1.0f * y,
                                        -0.5f + 0.3f * x + 0.7f * y);
  };

  size_t size_P = mesh.P.size();
  size_t size_T = mesh.T.size();
  Eigen::Matrix<float, Eigen::Dynamic, 1> u(2 * size_P + 2 * size_T);
  for (size_t i = 0; i < size_P; ++i) {
    femib::types::dvec<float, 2> uv = u_exact(mesh.P[i](0), mesh.P[i](1));
    u(static_cast<Eigen::Index>(i)) = uv(0);
    u(static_cast<Eigen::Index>(size_P + i)) = uv(1);
  }
  for (int n = 0; n < size_T; ++n) {
    femib::types::dvec<float, 2> c =
        femib::finite_element::find_center_of<float, 2>(mesh[n]);
    femib::types::dvec<float, 2> uv = u_exact(c(0), c(1));
    u(static_cast<Eigen::Index>(2 * size_P + n)) = uv(0);
    u(static_cast<Eigen::Index>(2 * size_P + size_T + n)) = uv(1);
  }

  std::vector<femib::types::dvec<float, 2>> query_points = {
      {0.31f, 0.42f}, {0.55f, 0.11f}, {0.72f, 0.83f}, {0.05f, 0.95f}};
  std::vector<femib::types::dvec<float, 2>> interpolated =
      femib::ib::interpolate_velocity<float, 2>(V, u, query_points);

  REQUIRE_EQ(interpolated.size(), query_points.size());
  for (size_t k = 0; k < query_points.size(); ++k) {
    femib::types::dvec<float, 2> expected =
        u_exact(query_points[k](0), query_points[k](1));
    CHECK_LT((interpolated[k] - expected).norm(), 1e-4f);
  }
}

TEST_CASE("spread_force conserves total force") {
  femib::types::mesh<float, 2> mesh = make_unit_square_mesh(6);
  mesh.init();

  femib::finite_element::finite_element<float, 2, 2> fe =
      femib::finite_element::create_finite_element_P1_B_2d2d<float, 2, 2>();
  femib::finite_element_space::finite_element_space V = {.finite_element = fe,
                                                         .mesh = mesh};
  V.nodes = fe.build_nodes(mesh);

  std::vector<femib::types::dvec<float, 2>> points = {
      {0.31f, 0.42f}, {0.55f, 0.11f}, {0.72f, 0.83f}};
  std::vector<femib::types::dvec<float, 2>> forces = {
      {1.0f, -0.5f}, {0.2f, 0.3f}, {-0.7f, 0.1f}};
  float dS = 0.05f;

  Eigen::Matrix<float, Eigen::Dynamic, 1> rhs =
      femib::ib::spread_force<float, 2>(V, points, forces, dS);

  size_t size_P = mesh.P.size();
  size_t size_T = mesh.T.size();
  REQUIRE_EQ(rhs.rows(), 2 * size_P + 2 * size_T);

  float total_x = 0.0f, total_y = 0.0f;
  for (size_t i = 0; i < size_P; ++i) {
    total_x += rhs(static_cast<Eigen::Index>(i));
    total_y += rhs(static_cast<Eigen::Index>(size_P + i));
  }
  for (size_t n = 0; n < size_T; ++n) {
    total_x += rhs(static_cast<Eigen::Index>(2 * size_P + n));
    total_y += rhs(static_cast<Eigen::Index>(2 * size_P + size_T + n));
  }

  float expected_x = 0.0f, expected_y = 0.0f;
  for (auto &force : forces) {
    expected_x += force(0) * dS;
    expected_y += force(1) * dS;
  }

  CHECK_LT(std::abs(total_x - expected_x), 1e-3f);
  CHECK_LT(std::abs(total_y - expected_y), 1e-3f);
}

TEST_CASE("testing force spreading of a given elastic ring") {

  femib::types::mesh<float, 2> mesh = make_unit_square_mesh(16);
  mesh.init();

  femib::finite_element::finite_element<float, 2, 2> fe =
      femib::finite_element::create_finite_element_P1_B_2d2d<float, 2, 2>();
  femib::finite_element_space::finite_element_space V = {.finite_element = fe,
                                                         .mesh = mesh};
  V.nodes = fe.build_nodes(mesh);

  femib::types::dvec<float, 2> center(0.5f, 0.5f);
  femib::ib::ring<float, 2> r =
      femib::ib::build_ring<float, 2>(center, 0.2f, 8, 3.0f);

  // Perturb away from the rest configuration so the force is non-trivial.
  for (size_t k = 0; k < r.X.size(); ++k) {
    r.X[k](0) += 0.01f * std::sin(3.0f * static_cast<float>(k));
    r.X[k](1) += 0.01f * std::cos(2.0f * static_cast<float>(k));
  }
  std::vector<femib::types::dvec<float, 2>> forces =
      femib::ib::elastic_force(r);

  Eigen::Matrix<float, Eigen::Dynamic, 1> spread_forces =
      femib::ib::spread_force<float, 2>(r, V);

  size_t size_P = mesh.P.size();
  size_t size_T = mesh.T.size();
  REQUIRE_EQ(spread_forces.rows(), 2 * size_P + 2 * size_T);

  float total_x = 0.0f, total_y = 0.0f;
  for (size_t i = 0; i < size_P; ++i) {
    total_x += spread_forces(static_cast<Eigen::Index>(i));
    total_y += spread_forces(static_cast<Eigen::Index>(size_P + i));
  }
  for (int n = 0; n < size_T; ++n) {
    total_x += spread_forces(static_cast<Eigen::Index>(2 * size_P + n));
    total_y +=
        spread_forces(static_cast<Eigen::Index>(2 * size_P + size_T + n));
  }

  float expected_x = 0.0f, expected_y = 0.0f;
  for (auto &force : forces) {
    expected_x += force(0) * r.dS;
    expected_y += force(1) * r.dS;
  }

  CHECK_LT(std::abs(total_x - expected_x), 1e-3f);
  CHECK_LT(std::abs(total_y - expected_y), 1e-3f);
}
