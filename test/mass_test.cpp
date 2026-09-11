#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../src/affine/affine.hpp"
#include "../src/finite_element/P1_2d1d.hpp"
#include "../src/finite_element/P1_2d2d.hpp"
#include "../src/finite_element_space/finite_element_space.hpp"
#include "../src/gauss/gauss.hpp"
#include "../src/gauss/gauss_lagrange_2_2d.hpp"
#include "../src/mesh/mesh.hpp"
#include <doctest/doctest.h>

#include "stokes_t.hpp"

TEST_CASE("testing mass_matrix") {
  femib::types::mesh<float, 2> mesh;
  mesh.P = {femib::types::dvec<float, 2>(0.0f, 0.0f),
            femib::types::dvec<float, 2>(1.0f, 0.0f),
            femib::types::dvec<float, 2>(0.0f, 1.0f)};
  mesh.T = {femib::types::ditrian<2>(0, 1, 2)};
  mesh.E = {0, 1, 2};
  mesh.init();

  femib::finite_element::finite_element<float, 2, 2> fe =
      femib::finite_element::create_finite_element_P1_2d2d<float, 2, 2>();
  femib::finite_element_space::finite_element_space V = {.finite_element = fe,
                                                         .mesh = mesh};
  V.nodes = fe.build_nodes(mesh);

  femib::gauss::rule<float, 2> rule =
      femib::gauss::create_gauss_2_2d<float, 2>();

  Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> M =
      femib::stokes_t::mass_matrix<float, 2>(V, rule);

  REQUIRE_EQ(M.rows(), 6);
  REQUIRE_EQ(M.cols(), 6);

  float diag = 1.0f / 12.0f;
  float off_same_component = 1.0f / 24.0f;

  // Same-component pairs (both x, indices 0-2, or both y, indices 3-5).
  int same_component_pairs[6][2] = {{0, 0}, {1, 1}, {2, 2},
                                    {3, 3}, {4, 4}, {5, 5}};
  for (auto &pr : same_component_pairs) {
    CHECK_LT(std::abs(M(pr[0], pr[1]) - diag), 1e-5f);
  }
  int off_diag_same_component[6][2] = {{0, 1}, {0, 2}, {1, 2},
                                       {3, 4}, {3, 5}, {4, 5}};
  for (auto &pr : off_diag_same_component) {
    CHECK_LT(std::abs(M(pr[0], pr[1]) - off_same_component), 1e-5f);
    CHECK_LT(std::abs(M(pr[1], pr[0]) - off_same_component),
             1e-5f); // symmetric
  }
  // Cross-component pairs (one x-basis-function, one y-basis-function) are
  // exactly zero
  for (int i = 0; i < 3; ++i) {
    for (int j = 3; j < 6; ++j) {
      CHECK_LT(std::abs(M(i, j)), 1e-6f);
      CHECK_LT(std::abs(M(j, i)), 1e-6f);
    }
  }
}
