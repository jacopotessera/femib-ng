#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../src/gauss/gauss.hpp"
#include "../src/mesh/mesh.hpp"
#include "../src/types/types.hpp"
#include <doctest/doctest.h>

const float EPSILON = std::numeric_limits<float>::epsilon();

TEST_CASE("testing mesh") {

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

  {
    std::function<float(femib::types::dvec<float, 2>)> f =
        [](femib::types::dvec<float, 2> x) { return 1; };
    float area = femib::mesh::integrate<float, 2>(rule, f, mesh);
    CHECK(std::fabs(area - 1.0) <= EPSILON);
  }
  {
    std::function<float(femib::types::dvec<float, 2>)> f =
        [](femib::types::dvec<float, 2> x) { return x(0) + x(1); };
    float integral = femib::mesh::integrate<float, 2>(rule, f, mesh);
    CHECK(std::fabs(integral - 1.0) <= EPSILON);
  }
  {
    std::function<float(femib::types::dvec<float, 2>)> f =
        [](femib::types::dvec<float, 2> x) { return x(0) + x(1) + 1; };
    float integral = femib::mesh::integrate<float, 2>(rule, f, mesh);
    CHECK(std::fabs(integral - 2.0) <= EPSILON);
  }

  femib::types::mesh<float, 2> mesh_from_file = femib::mesh::read<float, 2>(
      "../mesh/p0.mat", "../mesh/t0.mat", "../mesh/e0.mat");
  CHECK(std::fabs(mesh_from_file.P[0][0] - (-1.0)) <= EPSILON);
  CHECK(mesh_from_file.T[0][0] == 0);
  CHECK(mesh_from_file.P.size() == 5);
  CHECK(mesh_from_file.T.size() == 4);
  CHECK(mesh_from_file.E.size() == 4);
}
