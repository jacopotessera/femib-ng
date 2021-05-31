#ifndef FEMIB_HPP_INCLUDED_
#define FEMIB_HPP_INCLUDED_

#include "../affine/affine.hpp"
#include "../finite_element_space/finite_element_space.hpp"
#include "../gauss/gauss.hpp"
#include "../mesh/mesh.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <functional>
#include <vector>

int get_index(const femib::types::nodes<float, 2> &nodes, int i, int n) {
  return nodes.T[n][i];
}

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
std::function<T(femib::types::dvec<T, d>)>
ddot(femib::types::F<float, 2, 1> a, femib::types::F<float, 2, 1> b) {
  return [a, b](femib::types::dvec<float, 2> x) {
    return a.dx(x)[0] * b.dx(x)[0] + a.dx(x)[1] * b.dx(x)[1];
  };
}

template <typename T, int d, int e>
std::function<T(femib::types::dvec<T, d>)>
forz(femib::types::F<float, 2, 1> a) {
  return [a](femib::types::dvec<float, 2> x) {
    return -400 * x(0) * x(1) * a.x(x)[0];
  };
}

template <typename T, int d, int e>
std::vector<Eigen::Triplet<T>>
build_edges(femib::finite_element_space::finite_element_space<T, d, e> s,
            std::function<T(femib::types::dvec<T, d>)> b) {
  std::vector<Eigen::Triplet<float>> B;
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
// std::vector<Eigen::Triplet<T>>
MandF<T>
build_diagonal(femib::finite_element_space::finite_element_space<T, d, e> s,
               femib::gauss::rule<T, d> rule,
               std::function<std::function<T(femib::types::dvec<T, d>)>(
                   femib::types::F<float, 2, 1>, femib::types::F<float, 2, 1>)>
                   fff,
               std::function<std::function<T(femib::types::dvec<T, d>)>(
                   femib::types::F<float, 2, 1>)>
                   ggg) {
  std::vector<Eigen::Triplet<float>> MM;
  std::vector<Eigen::Triplet<float>> FF;
  for (int n = 0; n < s.mesh.T.size(); ++n) {
    femib::types::dtrian<float, 2> t = s.mesh[n];
    for (int i = 0; i < s.finite_element.base_functions.size(); ++i) {
      femib::types::F<float, 2, 1> a;

      a.x =
          [&](const femib::types::dvec<float, 2> &x) {
            return (s.finite_element.base_functions[i].x(
                femib::affine::affineBinv(t) *
                (x - femib::affine::affineb(t))));
          },
      a.dx = [&](const femib::types::dvec<float, 2> &x) {
        return (femib::affine::affineBinv(t) *
                s.finite_element.base_functions[i].dx(
                    femib::affine::affineBinv(t) *
                    (x - femib::affine::affineb(t))));
      };
      for (int j = 0; j < s.finite_element.base_functions.size(); ++j) {
        femib::types::F<float, 2, 1> b;
        b.dx = [&](const femib::types::dvec<float, 2> &x) {
          return (femib::affine::affineBinv(t) *
                  s.finite_element.base_functions[j].dx(
                      femib::affine::affineBinv(t) *
                      (x - femib::affine::affineb(t))));
        };
        float m = femib::mesh::integrate<float, 2>(rule, fff(a, b), t);
        MM.push_back(Eigen::Triplet<float>(get_index(s.nodes, i, n),
                                           get_index(s.nodes, j, n), m));
      }
      float f_ = femib::mesh::integrate<float, 2>(rule, ggg(a), t);
      FF.push_back(Eigen::Triplet<float>(get_index(s.nodes, i, n), 0, f_));
    }
  }
  return {MM, FF};
}

} // namespace femib::poisson
#endif
