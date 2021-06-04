#include "../src/affine/affine.hpp"
#include "../src/femib/poisson.hpp"
#include "../src/finite_element/P1_2d1d.hpp"
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
      mesh_dir + "p5.mat", mesh_dir + "t5.mat", mesh_dir + "e5.mat");
  mesh.init();

  femib::finite_element::finite_element<float, 2, 1> f =
      femib::finite_element::create_finite_element_P1_2d1d<float, 2, 1>();

  femib::finite_element_space::finite_element_space<float, 2, 1> s = {f, mesh};
  s.nodes = f.build_nodes(mesh);

  femib::poisson::poisson<float, 2, 1> poisson = {s};

  femib::poisson::init<float, 2, 1>(poisson, rule);

  Eigen::Matrix<float, Eigen::Dynamic, 1> xx =
      femib::poisson::solve<float, 2, 1>(poisson);

  std::for_each(s.nodes.T.begin(), s.nodes.T.end(),
                femib::poisson::print_node_generator<float, 2, 1>(s, xx));
}
