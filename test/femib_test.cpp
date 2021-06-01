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

  femib::finite_element::finite_element<float, 2, 1> f =
      femib::finite_element::create_finite_element_P1_2d1d<float, 2, 1>();

  femib::finite_element_space::finite_element_space<float, 2, 1> s = {f, mesh};
  s.nodes = f.build_nodes(mesh);

  femib::poisson::poisson<float, 2, 1> poisson = {s};

  /*
    femib::poisson::MandF<float> mandF =
        femib::poisson::build_diagonal<float, 2, 1>(
            s, rule, femib::poisson::ddot<float, 2, 1>,
            femib::poisson::forz<float, 2, 1>);

    std::vector<Eigen::Triplet<float>> M = mandF.M; // poisson.M(s, s);
    std::vector<Eigen::Triplet<float>> F = mandF.F; // poisson.f(s);

    std::function<float(femib::types::dvec<float, 2>)> b =
        [](const femib::types::dvec<float, 2> &x) {
          return 10 * x(0) * x(1) * x(1);
        };

    std::vector<Eigen::Triplet<float>> B =
        femib::poisson::build_edges<float, 2, 1>(s, b);*/

  std::vector<int> not_edges = femib::poisson::build_not_edges<float, 2, 1>(s);
  /*
      Eigen::Matrix<float, Eigen::Dynamic, 1> dB =
          femib::poisson::triplets2dense<float>(B, s.nodes.P.size(), 1);
      Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> dM =
          femib::poisson::triplets2dense<float>(M, s.nodes.P.size(),
                                                s.nodes.P.size());
      Eigen::Matrix<float, Eigen::Dynamic, 1> dF =
          femib::poisson::triplets2dense<float>(F, s.nodes.P.size(), 1);
    */

  femib::poisson::init<float, 2, 1>(poisson, rule);

  femib::poisson::AAAbbb<float> aaabbb = femib::poisson::remove_edges<float>(
      poisson.dM, poisson.dF, poisson.dB, s.nodes.P.size(), not_edges);

  Eigen::Matrix<float, Eigen::Dynamic, 1> xxx =
      aaabbb.AAA.colPivHouseholderQr().solve(aaabbb.bbb);

  Eigen::Matrix<float, Eigen::Dynamic, 1> xx = femib::poisson::add_edges<float>(
      xxx, poisson.dB, s.nodes.P.size(), not_edges, s.nodes.E);

  std::for_each(s.nodes.T.begin(), s.nodes.T.end(),
                femib::poisson::print_node_generator<float, 2, 1>(s, xx));
}
