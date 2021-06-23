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
#include <ctime>
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

int main() {

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

  // STOKES
  femib::stokes_t::stokes<float, 2> stokes;
  stokes.V = v;
  stokes.Q = q;
  femib::stokes_t::init<float, 2>(stokes, rule);

  femib::types::box<float, 2> box = femib::mesh::find_box<float, 2>(mesh);

  femib::types::box<float, 2> boxx =
      femib::mesh::lin_spaced<float, 2>(box, 0.1);

  std::string id = getTime();

  std::string dbname = "femib_test";
  femib::mongo::save_sim(dbname, id);

  // std::string id = "666";
  int TMAX = 100;
  for (int t = 0; t < TMAX; t++) {
    femib::stokes_t::advance<float, 2>(stokes);
    femib::mongo::plot_data p = {id, t, stokes.plot[t], {}, {}};
    femib::mongo::save_plot_data(dbname, p);
  }

  // Eigen::Matrix<float, Eigen::Dynamic, 1> xx = stokes.solution[TMAX - 1];

  // std::cout << xx << std::endl;
  // std::for_each(v.nodes.T.begin(), v.nodes.T.end(),
  //              femib::util::print_node_generator<float, 2, 2>(v, xx));
}
