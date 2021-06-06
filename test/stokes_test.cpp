#include "../src/affine/affine.hpp"
#include "../src/femib/stokes.hpp"
#include "../src/finite_element/P0_2d1d.hpp"
#include "../src/finite_element/P1+B_2d2d.hpp"
#include "../src/finite_element_space/finite_element_space.hpp"
#include "../src/gauss/gauss.hpp"
#include "../src/gauss/gauss_lagrange_2_2d.hpp"
#include "../src/mesh/mesh.hpp"
#include "../src/read/read.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <iostream>
#include <vector>

int main() {

  femib::gauss::rule<float, 2> rule =
      femib::gauss::create_gauss_2_2d<float, 2>();
  std::string mesh_dir = MESH_DIR;
  femib::types::mesh<float, 2> mesh = femib::mesh::read<float, 2>(
      mesh_dir + "p0.mat", mesh_dir + "t0.mat", mesh_dir + "e0.mat");
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
  femib::stokes::stokes<float, 2> stokes;
  stokes.V = v;
  stokes.Q = q;
  femib::stokes::init<float, 2>(stokes, rule);

  Eigen::Matrix<float, Eigen::Dynamic, 1> xx =
      femib::stokes::solve<float, 2, 1>(stokes);
  //std::cout << xx << std::endl;
  std::for_each(
      v.nodes.T.begin(), v.nodes.T.end(),
      femib::util::print_node_generator<float, 2, 2>(
          v, xx.block(0, 0, (v.nodes.P.size() - 2 * v.nodes.T.size()) / 2, 1)));

  std::for_each(
      v.nodes.T.begin(), v.nodes.T.end(),
      femib::util::print_node_generator<float, 2, 2>(
          v, xx.block((v.nodes.P.size() - 2 * v.nodes.T.size()) / 2, 0,
                      (v.nodes.P.size() - 2 * v.nodes.T.size()) / 2, 1)));
}
