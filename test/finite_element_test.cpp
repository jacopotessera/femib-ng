#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../src/finite_element/P1_2d1d.hpp"
#include "../src/finite_element/finite_element.hpp"
#include "../src/types/types.hpp"
#include <doctest/doctest.h>

TEST_CASE("testing finite_element") {

  femib::finite_element::finite_element<float, 2, 1> f =
      femib::finite_element::create_finite_element_P1_2d1d<float, 2, 1>();
  bool e = true;
  for (int i = 0; e && i < f.base_nodes.size(); i++) {
    for (int j = 0; e && j < f.base_functions.size(); j++) {
      e = e && (f.base_functions[j](f.base_nodes[i]).x[0] == (i == j ? 1 : 0));
    }
  }

  femib::types::mesh<float, 2> mesh = {
      {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}, {0.5, 0.5}},
      {{0, 1, 4}, {1, 2, 4}, {2, 3, 4}, {3, 0, 4}},
      {0, 1, 2, 3}};
  mesh.init();
  f.build_nodes(mesh);
  CHECK(e);
}
