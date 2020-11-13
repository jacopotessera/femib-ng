#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../src/affine/affine.hpp"
#include "../src/finite_element_space/finite_element_space.hpp"
#include "../src/gauss/gauss.hpp"
#include "../src/mesh/mesh.hpp"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <doctest/doctest.h>
#include <spdlog/spdlog.h>

int get_index(const femib::types::nodes<float, 2> &nodes, int i, int n) {
  return nodes.T[n][i];
}

TEST_CASE("testing finite_element_space") {

  std::vector<Eigen::Triplet<float>> M;
  std::vector<Eigen::Triplet<float>> F;

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
            [&](femib::types::dvec<float, 2> x) {
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
            [&](femib::types::dvec<float, 2> x) {
              float a_0 =
                  -1*(f.base_functions[i].x(affineBinv(t) * (x - affineb(t))))[0];
              return a_0;
            },
            t);
        //std::cout << "n: " << n << ", i: " << i << ", j: " << j
        //          << " -> m: " << m << std::endl;
        //std::cout << "n: " << n << ", i: " << i << ", j: " << j
        //          << " -> f: " << f_ << std::endl;
        M.push_back(Eigen::Triplet<float>(get_index(s.nodes, i, n),
                                          get_index(s.nodes, j, n), m));
        F.push_back(Eigen::Triplet<float>(get_index(s.nodes, i, n), 0, f_));
      }
    }
  }
  std::vector<Eigen::Triplet<float>> B;
  
  std::function<float(femib::types::dvec<float,2>)> b = [](const femib::types::dvec<float,2> &x) {
    return x(0)+x(1);
  };
  for(int e : s.nodes.E) {
    B.push_back(Eigen::Triplet<float>(e,0,b(s.nodes.P[e])));
  }

  std::vector<int> not_edges ;
  for(int i = 0; i< s.nodes.P.size(); i++){
    if(std::find(s.nodes.E.begin(),s.nodes.E.end(),i) == s.nodes.E.end()){
      not_edges.push_back(i);
    }
  }

  Eigen::SparseMatrix<float> sB = Eigen::SparseMatrix<float>(s.nodes.P.size(), 1);
  sB.setFromTriplets(B.begin(), B.end());
  Eigen::Matrix<float, Eigen::Dynamic,1> dB = Eigen::Matrix<float,Eigen::Dynamic,1>(sB);
  Eigen::SparseMatrix<float> sM = Eigen::SparseMatrix<float>(s.nodes.P.size(), s.nodes.P.size());
  sM.setFromTriplets(M.begin(), M.end());
  Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic > dM = Eigen::Matrix<float,Eigen::Dynamic, Eigen::Dynamic>(sM);
  Eigen::Matrix<float, Eigen::Dynamic,1> ss =  dM(not_edges,s.nodes.E)*dB(s.nodes.E,Eigen::all);
  Eigen::Matrix<float,Eigen::Dynamic,1> mf;
  mf.resize(s.nodes.P.size(),1);
  for(int i = 0; i< s.nodes.P.size(); i++){
    auto k = std::find(not_edges.begin(),not_edges.end(),i);
    mf(i,0) = 0.0;
    if(k != not_edges.end()){
      mf(i,0) = ss(k-not_edges.begin(),0);
    }
  }

  Eigen::SparseMatrix<float> sF = Eigen::SparseMatrix<float>(5, 1);
  sF.setFromTriplets(F.begin(), F.end());
  Eigen::Matrix<float, Eigen::Dynamic,1> dF = Eigen::Matrix<float,Eigen::Dynamic,1>(sF);
  Eigen::Matrix<float,Eigen::Dynamic,1> bbb = (dF - mf)(not_edges,0);
  Eigen::Matrix<float,Eigen::Dynamic,1> AAA = (dM)(not_edges,not_edges);
  std::cout << AAA.colPivHouseholderQr().solve(bbb) << std::endl;
  

  CHECK(true);
}
