#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include "../src/finite_element_space/finite_element_space.hpp"
#include "../src/gauss/gauss.hpp"

TEST_CASE("testing finite_element_space") {

  femib::gauss::node<float, 2> node1 = {1.0 / 6.0, {1.0 / 6.0, 1.0 / 6.0}};
  femib::gauss::node<float, 2> node2 = {1.0 / 6.0, {1.0 / 6.0, 2.0 / 3.0}};
  femib::gauss::node<float, 2> node3 = {1.0 / 6.0, {2.0 / 3.0, 1.0 / 6.0}};
  femib::gauss::rule<float, 2> rule = {{node1, node2, node3}};

  femib::types::mesh<float, 2> mesh = {
      // P
      {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}, {0.5, 0.5}},
      // T
      {{0, 1, 4}, {1, 2, 4}, {2, 3, 4}, {3, 0, 4}

      },
      // E
      {0, 1, 2, 3}};

  femib::finite_element::finite_element<float, 2, 1> f =
      femib::finite_element::create_finite_element_P1_2d1d<float, 2, 1>();



    CHECK(true);
}
