#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../src/affine/affine.hpp"
#include "../src/femib/stokes_t.hpp"
#include "../src/finite_element/P0_2d1d.hpp"
#include "../src/finite_element/P1+B_2d2d.hpp"
#include "../src/finite_element_space/finite_element_space.hpp"
#include "../src/gauss/gauss.hpp"
#include "../src/gauss/gauss_lagrange_2_2d.hpp"
#include "../src/mesh/mesh.hpp"
#include "../src/mongo/mongo.hpp"
#include "../src/read/read.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <doctest/doctest.h>
#include <iostream>
#include <stdio.h>
#include <sys/time.h>
#include <vector>

std::string getTime() {
  timeval curTime;

  gettimeofday(&curTime, NULL);

  int milli = curTime.tv_usec / 1000;
  char buf[sizeof "2011-10-08T07:07:09.000Z"];
  strftime(buf, sizeof buf, "%FT%T", gmtime(&curTime.tv_sec));
  sprintf(buf, "%s.%dZ", buf, milli);

  return buf;
}

namespace {

// Builds the mesh + V (velocity, P1+bubble) + Q (pressure, P0) + stokes_t
// struct, already through femib::stokes_t::init, shared by every TEST_CASE
// below -- factored out to avoid re-deriving the same boilerplate three
// times (mirrors the intent of test/stokes_test.cpp's own
// make_unit_square_mesh helper, Task 11).
femib::stokes_t::stokes<float, 2>
make_stokes_t_fixture(femib::types::mesh<float, 2> &mesh_out) {
  femib::gauss::rule<float, 2> rule =
      femib::gauss::create_gauss_2_2d<float, 2>();
  std::string mesh_dir = MESH_DIR;
  femib::types::mesh<float, 2> mesh = femib::mesh::read<float, 2>(
      mesh_dir + "p3.mat", mesh_dir + "t3.mat", mesh_dir + "e3.mat");
  mesh.init();

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

  femib::stokes_t::stokes<float, 2> s;
  s.V = v;
  s.Q = q;
  femib::stokes_t::init<float, 2>(s, rule);

  mesh_out = mesh;
  return s;
}

} // namespace

TEST_CASE("testing femib stokes_t pipeline with mongo persistence") {
  femib::types::mesh<float, 2> mesh;
  femib::stokes_t::stokes<float, 2> stokes = make_stokes_t_fixture(mesh);

  femib::types::box<float, 2> box = femib::mesh::find_box<float, 2>(mesh);

  femib::types::box<float, 2> boxx =
      femib::mesh::lin_spaced<float, 2>(box, 0.1);

  std::string id = getTime();

  std::string dbname = "femib_test";
  femib::mongo::save_sim(dbname, id);

  int TMAX = 100;
  for (int t = 0; t < TMAX; t++) {
    femib::stokes_t::advance<float, 2>(stokes);
    femib::mongo::plot_data p = {id, t, stokes.plot[t], {}, {}};
    femib::mongo::save_plot_data(dbname, p);
  }

  CHECK(stokes.solution.size() == (size_t)TMAX);
  CHECK(stokes.solution.back().allFinite());
}

TEST_CASE("testing advance") {
  femib::types::mesh<float, 2> mesh;
  femib::stokes_t::stokes<float, 2> s = make_stokes_t_fixture(mesh);

  femib::stokes_t::advance<float, 2>(s);
  CHECK_NOTHROW(femib::stokes_t::advance<float, 2>(s));

  CHECK(s.solution.back().rows() == s.V.nodes.P.size() + s.Q.nodes.P.size());
  CHECK(s.solution.back().topRows(s.V.nodes.P.size()).allFinite());
}

TEST_CASE("testing advance over several timesteps") {
  femib::types::mesh<float, 2> mesh;
  femib::stokes_t::stokes<float, 2> s = make_stokes_t_fixture(mesh);

  int rowsV = s.V.nodes.P.size();
  const int n_steps = 8;
  std::vector<float> velocity_norm;
  std::vector<float> step_diff_norm;

  Eigen::Matrix<float, Eigen::Dynamic, 1> prev_v;
  for (int step = 0; step < n_steps; ++step) {
    femib::stokes_t::advance<float, 2>(s);
    Eigen::Matrix<float, Eigen::Dynamic, 1> v =
        s.solution.back().topRows(rowsV);
    CHECK(v.allFinite());
    velocity_norm.push_back(v.norm());
    if (step > 0) {
      step_diff_norm.push_back((v - prev_v).norm());
    }
    prev_v = v;
  }

  float max_norm =
      *std::max_element(velocity_norm.begin(), velocity_norm.end());
  float min_norm =
      *std::min_element(velocity_norm.begin(), velocity_norm.end());
  CHECK(max_norm < 10.0f * min_norm);

  CHECK(step_diff_norm.back() < 1.5f * step_diff_norm.front());

  Eigen::Matrix<float, Eigen::Dynamic, 1> p =
      s.solution.back().bottomRows(s.Q.nodes.P.size());
  CHECK(std::abs(s.domain_integral_row.dot(p)) < 1e-4f);
}
