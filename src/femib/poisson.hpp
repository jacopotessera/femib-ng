#ifndef FEMIB_POISSON_HPP_INCLUDED_
#define FEMIB_POISSON_HPP_INCLUDED_

#include "../affine/affine.hpp"
#include "../finite_element_space/finite_element_space.hpp"
#include "../gauss/gauss.hpp"
#include "../mesh/mesh.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <functional>
#include <iostream>
#include <vector>

namespace femib::poisson {

template <typename T, int d, int e> struct poisson {

  femib::finite_element_space::finite_element_space<T, d, e> V;
  std::function<std::vector<Eigen::Triplet<T>>(
      femib::finite_element_space::finite_element_space<T, d, e> a,
      femib::finite_element_space::finite_element_space<T, d, e> b)>
      M;
  std::function<std::vector<Eigen::Triplet<T>>(
      femib::finite_element_space::finite_element_space<T, d, e> a)>
      f;

  Eigen::Matrix<T, Eigen::Dynamic, 1> dB;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> dM;
  Eigen::Matrix<T, Eigen::Dynamic, 1> dF;
};

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
triplets2dense(std::vector<Eigen::Triplet<T>> M, int rows, int cols) {

  Eigen::SparseMatrix<T> sM = Eigen::SparseMatrix<T>(rows, cols);
  sM.setFromTriplets(M.begin(), M.end());
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> dM =
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(sM);
  return dM;
}

template <typename T, int d, int e>
std::function<T(femib::types::dvec<T, d>)> ddot(femib::types::F<T, d, e> a,
                                                femib::types::F<T, d, e> b) {
  return [a, b](femib::types::dvec<T, d> x) {
    return a.dx(x)[0] * b.dx(x)[0] + a.dx(x)[1] * b.dx(x)[1];
  };
}

template <typename T, int d, int e>
std::function<T(femib::types::dvec<T, d>)> forz(femib::types::F<T, d, e> a) {
  return [a](femib::types::dvec<T, d> x) {
    return -400 * x(0) * x(1) * a.x(x)[0];
  };
}

template <typename T, int d, int e>
std::vector<Eigen::Triplet<T>>
build_edges(femib::finite_element_space::finite_element_space<T, d, e> s,
            std::function<T(femib::types::dvec<T, d>)> b) {
  std::vector<Eigen::Triplet<T>> B;
  for (int i : s.nodes.E) {
    B.push_back(Eigen::Triplet<T>(i, 0, b(s.nodes.P[i])));
  }
  return B;
}

template <typename T, int d, int e>
std::vector<int>
build_not_edges(femib::finite_element_space::finite_element_space<T, d, e> s) {
  std::vector<int> not_edges;
  for (int i = 0; i < s.nodes.P.size(); i++) {
    if (std::find(s.nodes.E.begin(), s.nodes.E.end(), i) == s.nodes.E.end()) {
      not_edges.push_back(i);
    }
  }
  return not_edges;
}

template <typename T> struct MandF {
  std::vector<Eigen::Triplet<T>> M;
  std::vector<Eigen::Triplet<T>> F;
};

template <typename T, int d, int e>
MandF<T>
build_diagonal(femib::finite_element_space::finite_element_space<T, d, e> s,
               femib::gauss::rule<T, d> rule,
               std::function<std::function<T(femib::types::dvec<T, d>)>(
                   femib::types::F<T, d, e>, femib::types::F<T, d, e>)>
                   fff,
               std::function<std::function<T(femib::types::dvec<T, d>)>(
                   femib::types::F<T, d, e>)>
                   ggg) {
  std::vector<Eigen::Triplet<T>> MM;
  std::vector<Eigen::Triplet<T>> FF;
  for (int n = 0; n < s.mesh.T.size(); ++n) {
    femib::types::dtrian<T, d> t = s.mesh[n];
    for (int i = 0; i < s.finite_element.base_functions.size(); ++i) {
      femib::types::F<T, d, e> a;

      a.x =
          [&](const femib::types::dvec<T, d> &x) {
            return (s.finite_element.base_functions[i].x(
                femib::affine::affineBinv(t) *
                (x - femib::affine::affineb(t))));
          },
      a.dx = [&](const femib::types::dvec<T, d> &x) {
        return (femib::affine::affineBinv(t) *
                s.finite_element.base_functions[i].dx(
                    femib::affine::affineBinv(t) *
                    (x - femib::affine::affineb(t))));
      };
      for (int j = 0; j < s.finite_element.base_functions.size(); ++j) {
        femib::types::F<T, d, e> b;
        b.dx = [&](const femib::types::dvec<T, d> &x) {
          return (femib::affine::affineBinv(t) *
                  s.finite_element.base_functions[j].dx(
                      femib::affine::affineBinv(t) *
                      (x - femib::affine::affineb(t))));
        };
        T m = femib::mesh::integrate<T, d>(rule, fff(a, b), t);
        MM.push_back(
            Eigen::Triplet<T>(femib::types::get_index<T, d>(s.nodes, i, n),
                              femib::types::get_index<T, d>(s.nodes, j, n), m));
      }
      T f_ = femib::mesh::integrate<T, d>(rule, ggg(a), t);
      FF.push_back(Eigen::Triplet<T>(
          femib::types::get_index<T, d>(s.nodes, i, n), 0, f_));
    }
  }
  return {MM, FF};
}

template <typename T> struct AAAbbb {
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> AAA;
  Eigen::Matrix<T, Eigen::Dynamic, 1> bbb;
};

template <typename T>
AAAbbb<T> remove_edges(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> dM,
                       Eigen::Matrix<T, Eigen::Dynamic, 1> dF,
                       Eigen::Matrix<T, Eigen::Dynamic, 1> dB, int rows,
                       std::vector<int> not_edges) {

  Eigen::Matrix<T, Eigen::Dynamic, 1> ss =
      dM(Eigen::all, Eigen::all) * dB(Eigen::all, Eigen::all);
  Eigen::Matrix<T, Eigen::Dynamic, 1> mf;
  mf.resize(rows, 1);
  for (int i = 0; i < rows; i++) {
    auto k = std::find(not_edges.begin(), not_edges.end(), i);
    mf(i, 0) = 0.0;
    if (k != not_edges.end()) {
      mf(i, 0) = ss(k - not_edges.begin(), 0);
    }
  }

  Eigen::Matrix<T, Eigen::Dynamic, 1> bbb = (dF - ss)(not_edges, 0);
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> AAA =
      (dM)(not_edges, not_edges);

  return {AAA, bbb};
}

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> add_edges(

    Eigen::Matrix<T, Eigen::Dynamic, 1> xxx,
    Eigen::Matrix<T, Eigen::Dynamic, 1> dB, int rows,
    std::vector<int> not_edges, std::vector<int> nodesE) {

  Eigen::Matrix<T, Eigen::Dynamic, 1> xx;
  xx.resize(rows, 1);

  for (int i = 0; i < rows; i++) {
    xx(i, 0) = 0.0;
    auto k = std::find(not_edges.begin(), not_edges.end(), i);
    if (k != not_edges.end()) {
      xx(i, 0) = xxx(k - not_edges.begin(), 0);
    }
    auto kk = std::find(nodesE.begin(), nodesE.end(), i);
    if (kk != nodesE.end()) {
      xx(i, 0) = dB(i, 0);
    }
  }

  return xx;
}

template <typename T, int d, int e>
auto print_node_generator(
    femib::finite_element_space::finite_element_space<T, d, e> s,
    Eigen::Matrix<T, Eigen::Dynamic, 1> xx) {

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

template <typename T, int d, int e>
void init(poisson<T, d, e> &s, femib::gauss::rule<T, d> &rule) {

  MandF<T> mandF =
      build_diagonal<T, d, e>(s.V, rule, ddot<T, d, e>, forz<T, d, e>);

  std::vector<Eigen::Triplet<T>> M = mandF.M; // poisson.M(s, s);
  std::vector<Eigen::Triplet<T>> F = mandF.F; // poisson.f(s);

  std::function<T(femib::types::dvec<T, d>)> b =
      [](const femib::types::dvec<T, d> &x) { return 10 * x(0) * x(1) * x(1); };

  std::vector<Eigen::Triplet<T>> B = build_edges<T, d, e>(s.V, b);

  std::vector<int> not_edges = build_not_edges<T, d, e>(s.V);

  s.dB = triplets2dense<T>(B, s.V.nodes.P.size(), 1);
  s.dM = triplets2dense<T>(M, s.V.nodes.P.size(), s.V.nodes.P.size());
  s.dF = triplets2dense<T>(F, s.V.nodes.P.size(), 1);
}

} // namespace femib::poisson
#endif
