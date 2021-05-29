//#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../src/affine/affine.hpp"
#include "../src/finite_element/P1_2d1d.hpp"
#include "../src/finite_element_space/finite_element_space.hpp"
#include "../src/gauss/gauss.hpp"
#include "../src/gauss/gauss_lagrange_2_2d.hpp"
#include "../src/mesh/mesh.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
//#include <doctest/doctest.h>
#include <iostream>
#include <vector>

using namespace femib::affine;

int get_index(const femib::types::nodes<float, 2> &nodes, int i, int n) {
  return nodes.T[n][i];
}

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

// TEST_CASE("testing finite_element_space") {
int main() {

  std::vector<Eigen::Triplet<float>> M;
  std::vector<Eigen::Triplet<float>> F;

  femib::gauss::rule<float, 2> rule =
      femib::gauss::create_gauss_2_2d<float, 2>();

  std::string mesh_dir = MESH_DIR;
  femib::types::mesh<float, 2> mesh = femib::mesh::read<float, 2>(
      mesh_dir + "p5.mat", mesh_dir + "t5.mat", mesh_dir + "e5.mat");

  // for(femib::types::dtrian<float,2> v : mesh){
  //  std::cout << v[0] << std::endl;
  //  std::cout << v[1] << std::endl;
  //  std::cout << v[2] << std::endl;
  //  std::cout << std::endl;
  //}

  femib::finite_element::finite_element<float, 2, 1> f =
      femib::finite_element::create_finite_element_P1_2d1d<float, 2, 1>();

  femib::finite_element_space::finite_element_space<float, 2, 1> s = {f, mesh};
  s.nodes = f.build_nodes(mesh);

  for (int n = 0; n < s.mesh.T.size(); ++n) {
    for (int i = 0; i < s.finite_element.base_functions.size(); ++i) {
      for (int j = 0; j < s.finite_element.base_functions.size(); ++j) {
        femib::types::F<float, 2, 1> a;
        femib::types::F<float, 2, 1> b;
        femib::types::dtrian<float, 2> t = s.mesh[n];
        a.dx = [&](const femib::types::dvec<float, 2> &x) {
          return (affineBinv(t) *
                  f.base_functions[i].dx(affineBinv(t) * (x - affineb(t))));
        };
        b.dx = [&](const femib::types::dvec<float, 2> &x) {
          return (affineBinv(t) *
                  f.base_functions[j].dx(affineBinv(t) * (x - affineb(t))));
        };
        // float m = femib::mesh::integrate<float,2>(rule,
        // [](femib::types::dvec<float,2>){return (float)1.0;}, t);
        float m = femib::mesh::integrate<float, 2>(
            rule,
            [&t, &f, i, j](femib::types::dvec<float, 2> x) {
              float a_0 =
                  (affineBinv(t) *
                   f.base_functions[i].dx(affineBinv(t) * (x - affineb(t))))[0];
              float a_1 =
                  (affineBinv(t) *
                   f.base_functions[i].dx(affineBinv(t) * (x - affineb(t))))[1];
              float b_0 =
                  (affineBinv(t) *
                   f.base_functions[j].dx(affineBinv(t) * (x - affineb(t))))[0];
              float b_1 =
                  (affineBinv(t) *
                   f.base_functions[j].dx(affineBinv(t) * (x - affineb(t))))[1];
              return a_0 * b_0 + a_1 * b_1;
              // return a.dx(x)[0]*b.dx(x)[0] + a.dx(x)[1]*b.dx(x)[1];
            },
            t);
        float f_ = femib::mesh::integrate<float, 2>(
            rule,
            [&t, &f, i](femib::types::dvec<float, 2> x) {
              float a_0 =
                  -40 * x(0) * x(1) *
                  (f.base_functions[i].x(affineBinv(t) * (x - affineb(t))))[0];
              return a_0;
            },
            t);
        // std::cout << "n: " << n << ", i: " << i << ", j: " << j
        //          << " -> m: " << m << std::endl;
        // std::cout << "n: " << n << ", i: " << i << ", j: " << j
        //          << " -> f: " << f_ << std::endl;
        M.push_back(Eigen::Triplet<float>(get_index(s.nodes, i, n),
                                          get_index(s.nodes, j, n), m));
        F.push_back(Eigen::Triplet<float>(get_index(s.nodes, i, n), 0, f_));
      }
    }
  }
  std::vector<Eigen::Triplet<float>> B;

  std::function<float(femib::types::dvec<float, 2>)> b =
      [](const femib::types::dvec<float, 2> &x) { return x(0) + x(1); };
  //[](const femib::types::dvec<float, 2> &x) { return 1.0; };
  for (int e : s.nodes.E) {
    // std::cout << e << ": " << b(s.nodes.P[e]) << std::endl;
    B.push_back(Eigen::Triplet<float>(e, 0, b(s.nodes.P[e])));
  }

  std::vector<int> not_edges;
  for (int i = 0; i < s.nodes.P.size(); i++) {
    if (std::find(s.nodes.E.begin(), s.nodes.E.end(), i) == s.nodes.E.end()) {
      not_edges.push_back(i);
    }
  }

  Eigen::SparseMatrix<float> sB =
      Eigen::SparseMatrix<float>(s.nodes.P.size(), 1);
  sB.setFromTriplets(B.begin(), B.end());
  Eigen::Matrix<float, Eigen::Dynamic, 1> dB =
      Eigen::Matrix<float, Eigen::Dynamic, 1>(sB);
  Eigen::SparseMatrix<float> sM =
      Eigen::SparseMatrix<float>(s.nodes.P.size(), s.nodes.P.size());
  sM.setFromTriplets(M.begin(), M.end());
  Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> dM =
      Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>(sM);
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

  Eigen::SparseMatrix<float> sF =
      Eigen::SparseMatrix<float>(s.nodes.P.size(), 1);
  sF.setFromTriplets(F.begin(), F.end());
  Eigen::Matrix<float, Eigen::Dynamic, 1> dF =
      Eigen::Matrix<float, Eigen::Dynamic, 1>(sF);
  Eigen::Matrix<float, Eigen::Dynamic, 1> bbb = (dF - ss)(not_edges, 0);
  Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> AAA =
      (dM)(not_edges, not_edges);
  Eigen::Matrix<float, Eigen::Dynamic, 1> xxx =
      AAA.colPivHouseholderQr().solve(bbb);
  Eigen::Matrix<float, Eigen::Dynamic, 1> xx;
  xx.resize(s.nodes.P.size(), 1);

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
  // CHECK(true);
}
