#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../src/read/read.hpp"
#include "name_reporter.h"
#include <cstdio>
#include <doctest/doctest.h>
#include <fstream>

TEST_CASE("testing read_mesh_file: file not found") {
  CHECK_THROWS(read_mesh_file<float, 2>("not-existing.mat", "not-existing.mat",
                                        "not-existing.mat"));
}

TEST_CASE("testing read_mesh_file") {
  std::string mesh_dir = MESH_DIR;
  const std::string p_path = mesh_dir + "p5.mat";
  const std::string t_path = mesh_dir + "t5.mat";
  const std::string e_path = mesh_dir + "e5.mat";

  femib::types::mesh<float, 2> mesh =
      read_mesh_file<float, 2>(p_path, t_path, e_path);

  CHECK(mesh.P[77](0) == doctest::Approx(-0.84436));
  CHECK(mesh.P[77](1) == doctest::Approx(0.84098));

  CHECK(mesh.T[1](0) == doctest::Approx(1));
  CHECK(mesh.T[1](1) == doctest::Approx(1921));
  CHECK(mesh.T[1](2) == doctest::Approx(692));

  CHECK(mesh.E[34] == doctest::Approx(63));
}
