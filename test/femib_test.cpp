#include "../src/affine/affine.hpp"
#include "../src/femib/femib.hpp"
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

auto printNode_generator(
    femib::finite_element_space::finite_element_space<float, 2, 1> s,
    Eigen::Matrix<float, Eigen::Dynamic, 1> xx) {

  return [&s, &xx](std::vector<int> t) {
    int j_0 = t[0];
    int j_1 = t[1];
    int j_2 = t[2];
    std::cout << s.nodes.P[j_0](0) << "\t" << s.nodes.P[j_0](1) << "\t"
              << xx(j_0) << std::endl;
    std::cout << s.nodes.P[j_1](0) << "\t" << s.nodes.P[j_1](1) << "\t"
              << xx(j_1) << std::endl;
    std::cout << s.nodes.P[j_2](0) << "\t" << s.nodes.P[j_2](1) << "\t"
              << xx(j_2) << std::endl;
    std::cout << s.nodes.P[j_0](0) << "\t" << s.nodes.P[j_0](1) << "\t"
              << xx(j_0) << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
  };
}

int main() {

  femib::gauss::rule<float, 2> rule =
      femib::gauss::create_gauss_2_2d<float, 2>();
  std::string mesh_dir = MESH_DIR;
  femib::types::mesh<float, 2> mesh = femib::mesh::read<float, 2>(
      mesh_dir + "p4.mat", mesh_dir + "t4.mat", mesh_dir + "e4.mat");

  femib::finite_element::finite_element<float, 2, 1> f =
      femib::finite_element::create_finite_element_P1_2d1d<float, 2, 1>();

  femib::finite_element_space::finite_element_space<float, 2, 1> s = {f, mesh};
  s.nodes = f.build_nodes(mesh);

  femib::poisson::poisson<float, 2, 1> poisson;
  poisson.V = s;

  femib::poisson::MandF<float> mandF =
      femib::poisson::build_diagonal<float, 2, 1>(
          s, rule, femib::poisson::ddot<float, 2, 1>,
          femib::poisson::forz<float, 2, 1>);

  std::vector<Eigen::Triplet<float>> M = mandF.M; // poisson.M(s, s);
  std::vector<Eigen::Triplet<float>> F = mandF.F; // poisson.f(s);

  std::function<float(femib::types::dvec<float, 2>)> b =
      [](const femib::types::dvec<float, 2> &x) { return 10 * x(0); };

  std::vector<Eigen::Triplet<float>> B =
      femib::poisson::build_edges<float, 2, 1>(s, b);

  std::vector<int> not_edges = femib::poisson::build_not_edges<float, 2, 1>(s);

  Eigen::Matrix<float, Eigen::Dynamic, 1> dB =
      femib::poisson::triplets2dense<float>(B, s.nodes.P.size(), 1);
  Eigen::Matrix<float, Eigen::Dynamic, 1> dM =
      femib::poisson::triplets2dense<float>(M, s.nodes.P.size(),
                                            s.nodes.P.size());
  Eigen::Matrix<float, Eigen::Dynamic, 1> dF =
      femib::poisson::triplets2dense<float>(F, s.nodes.P.size(), 1);

  // cosa fa questo ???
  Eigen::Matrix<float, Eigen::Dynamic, 1> ss =
      dM(Eigen::all, Eigen::all) * dB(Eigen::all, Eigen::all);
  Eigen::Matrix<float, Eigen::Dynamic, 1> mf;
  mf.resize(s.nodes.P.size(), 1);
  for (int i = 0; i < s.nodes.P.size(); i++) {
    auto k = std::find(not_edges.begin(), not_edges.end(), i);
    mf(i, 0) = 0.0;
    if (k != not_edges.end()) {
      mf(i, 0) = ss(k - not_edges.begin(), 0);
    }
  }

  // cosa fa questo ???
  Eigen::Matrix<float, Eigen::Dynamic, 1> bbb = (dF - ss)(not_edges, 0);
  Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> AAA =
      (dM)(not_edges, not_edges);
  Eigen::Matrix<float, Eigen::Dynamic, 1> xxx =
      AAA.colPivHouseholderQr().solve(bbb);
  Eigen::Matrix<float, Eigen::Dynamic, 1> xx;
  xx.resize(s.nodes.P.size(), 1);

  // cosa fa questo ???
  for (int i = 0; i < s.nodes.P.size(); i++) {
    xx(i, 0) = 0.0;
    auto k = std::find(not_edges.begin(), not_edges.end(), i);
    if (k != not_edges.end()) {
      xx(i, 0) = xxx(k - not_edges.begin(), 0);
    }
    auto kk = std::find(s.nodes.E.begin(), s.nodes.E.end(), i);
    if (kk != s.nodes.E.end()) {
      xx(i, 0) = dB(i, 0);
    }
  }

  if (true) {
    std::for_each(s.nodes.T.begin(), s.nodes.T.end(),
                  printNode_generator(s, xx));
  }
}
