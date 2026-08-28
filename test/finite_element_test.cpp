#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../src/finite_element/P0_2d1d.hpp"
#include "../src/finite_element/P1+B_2d2d.hpp"
#include "../src/finite_element/P1_2d1d.hpp"
#include "../src/finite_element/P1_2d2d.hpp"
#include "../src/finite_element/finite_element.hpp"
#include "../src/types/types.hpp"
#include <doctest/doctest.h>

TEST_CASE("testing finite_element P1_2d1d") {
  femib::finite_element::finite_element<float, 2, 1> f =
      femib::finite_element::create_finite_element_P1_2d1d<float, 2, 1>();
  for (int i = 0; i < f.base_nodes.size(); ++i) {
    for (int j = 0; j < f.base_functions.size(); ++j) {
      CAPTURE(i);
      CAPTURE(j);
      float value = f.base_functions[j].x(f.base_nodes[i])(0);
      CHECK(value == doctest::Approx(i == j ? 1.0 : 0.0));
    }
  }
}

TEST_CASE("testing P1_2d1d build_nodes") {
  femib::finite_element::finite_element<float, 2, 1> f =
      femib::finite_element::create_finite_element_P1_2d1d<float, 2, 1>();
  femib::types::mesh<float, 2> mesh = {
      .P = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}, {0.5, 0.5}},
      .T = {{0, 1, 4}, {1, 2, 4}, {2, 3, 4}, {3, 0, 4}},
      .E = {0, 1, 2, 3}};
  mesh.init();
  femib::types::nodes<float, 2> nodes = f.build_nodes(mesh);
  CHECK(nodes.P.size() == mesh.P.size());
  CHECK(nodes.T.size() == mesh.T.size());
}

TEST_CASE("testing finite_element P0_2d1d") {
  femib::finite_element::finite_element<float, 2, 1> f =
      femib::finite_element::create_finite_element_P0_2d1d<float, 2, 1>();
  for (int i = 0; i < f.base_nodes.size(); ++i) {
    for (int j = 0; j < f.base_functions.size(); ++j) {
      CAPTURE(i);
      CAPTURE(j);
      float value = f.base_functions[j].x(f.base_nodes[i])(0);
      CHECK(value == doctest::Approx(i == j ? 1.0 : 0.0));
    }
  }
}

TEST_CASE("testing finite_element P1_2d2d") {
  femib::finite_element::finite_element<float, 2, 2> f =
      femib::finite_element::create_finite_element_P1_2d2d<float, 2, 2>();
  int n_nodes = f.base_nodes.size();
  CHECK(n_nodes == 3);
  CHECK(f.base_functions.size() == 6);
  for (int k = 0; k < f.base_functions.size(); ++k) {
    int component = k / n_nodes;
    int node_idx = k % n_nodes;
    for (int i = 0; i < n_nodes; ++i) {
      femib::types::dvec<float, 2> value =
          f.base_functions[k].x(f.base_nodes[i]);
      for (int c = 0; c < 2; ++c) {
        float expected = (c == component && i == node_idx) ? 1.0f : 0.0f;
        CHECK(value(c) == doctest::Approx(expected));
      }
    }
  }
}

TEST_CASE("testing finite_element P1+B_2d2d") {
  femib::finite_element::finite_element<float, 2, 2> f =
      femib::finite_element::create_finite_element_P1_B_2d2d<float, 2, 2>();
  int n_nodes = f.base_nodes.size();
  CHECK(n_nodes == 4);
  CHECK(f.base_functions.size() == 8);
  int component[8] = {0, 0, 0, 1, 1, 1, 0, 1};
  int node_idx[8] = {0, 1, 2, 0, 1, 2, 3, 3};
  for (int k = 0; k < f.base_functions.size(); ++k) {
    for (int i = 0; i < n_nodes; ++i) {
      femib::types::dvec<float, 2> value =
          f.base_functions[k].x(f.base_nodes[i]);
      for (int c = 0; c < 2; ++c) {
        float expected = (c == component[k] && i == node_idx[k]) ? 1.0f : 0.0f;
        CHECK(value(c) == doctest::Approx(expected));
      }
    }
  }
}
