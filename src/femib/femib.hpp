#ifndef FEMIB_HPP_INCLUDED_
#define FEMIB_HPP_INCLUDED_

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>

#include "../affine/affine.hpp"
#include "../finite_element_space/finite_element_space.hpp"
#include "../gauss/gauss.hpp"
#include "../mesh/mesh.hpp"
#include "../types/differential_operation.hpp"
#include "../types/types.hpp"

namespace femib::util {

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
triplets2dense(std::vector<Eigen::Triplet<T>> triplets, int rows, int cols) {
  Eigen::SparseMatrix<T> sparse_matrix = Eigen::SparseMatrix<T>(rows, cols);
  sparse_matrix.setFromTriplets(triplets.begin(), triplets.end());
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> dense_matrix =
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(sparse_matrix);
  return dense_matrix;
}

template <typename T, int d, int e>
femib::types::F<T, d, e> base_function2real_function(
    const femib::finite_element_space::finite_element_space<T, d, e> &v, int n,
    int i) {
  femib::types::F<T, d, e> a;
  a.x =
      [&v, n, i](const femib::types::dvec<T, d> &x) {
        return (v.finite_element.base_functions[i].x(
            femib::affine::affineBinv(v.mesh[n]) *
            (x - femib::affine::affineb(v.mesh[n]))));
      },
  a.dx = [&v, n, i](const femib::types::dvec<T, d> &x) {
    return (femib::affine::affineBinv(v.mesh[n]) *
            v.finite_element.base_functions[i].dx(
                femib::affine::affineBinv(v.mesh[n]) *
                (x - femib::affine::affineb(v.mesh[n]))));
  };
  return a;
}

template <typename T> struct build_diagonal_result {
  std::vector<Eigen::Triplet<T>> M;
  std::vector<Eigen::Triplet<T>> F;
};

template <typename T> struct solvable_equations {
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A;
  Eigen::Matrix<T, Eigen::Dynamic, 1> b;
};

template <typename T, int d, int e>
build_diagonal_result<T> build_diagonal(
    const femib::finite_element_space::finite_element_space<T, d, e> &v,
    const femib::gauss::rule<T, d> &rule,
    const std::function<std::function<T(femib::types::dvec<T, d>)>(
        femib::types::F<T, d, e>, femib::types::F<T, d, e>)> &fff,
    const std::function<std::function<T(femib::types::dvec<T, d>)>(
        femib::types::F<T, d, e>)> &ggg) {
  std::vector<Eigen::Triplet<T>> BB;
  std::vector<Eigen::Triplet<T>> FF;
  for (int n = 0; n < v.mesh.T.size(); ++n) {
    femib::types::dtrian<T, d> t = v.mesh[n];
    for (int i = 0; i < v.finite_element.base_functions.size(); ++i) {
      femib::types::F<T, d, e> a =
          femib::util::base_function2real_function<T, d, e>(v, n, i);
      for (int j = 0; j < v.finite_element.base_functions.size(); ++j) {
        femib::types::F<T, d, e> b =
            femib::util::base_function2real_function<T, d, e>(v, n, j);
        T m = femib::mesh::integrate<T, d>(rule, fff(a, b), t);
        BB.push_back(Eigen::Triplet<T>(v.nodes.get_index(i, n),
                                       v.nodes.get_index(j, n), m));
      }
      T f_ = femib::mesh::integrate<T, d>(rule, ggg(a), t);
      FF.push_back(Eigen::Triplet<T>(v.nodes.get_index(i, n), 0, f_));
    }
  }
  return {BB, FF};
}

template <typename T, int d>
std::vector<Eigen::Triplet<T>> build_non_diagonal(
    const femib::finite_element_space::finite_element_space<T, d, d> &v,
    const femib::finite_element_space::finite_element_space<T, d, 1> &q,
    const femib::gauss::rule<T, d> &rule,
    const std::function<std::function<T(femib::types::dvec<T, d>)>(
        femib::types::F<T, d, d>, femib::types::F<T, d, 1>)> &fff) {
  std::vector<Eigen::Triplet<T>> BB;
  for (int n = 0; n < v.mesh.T.size(); ++n) {
    femib::types::dtrian<T, d> t = v.mesh[n];
    for (int i = 0; i < v.finite_element.base_functions.size(); ++i) {
      femib::types::F<T, d, d> a =
          femib::util::base_function2real_function<T, d, d>(v, n, i);
      for (int j = 0; j < q.finite_element.base_functions.size(); ++j) {
        femib::types::F<T, d, 1> b =
            femib::util::base_function2real_function<T, d, 1>(q, n, j);
        T m = femib::mesh::integrate<T, d>(rule, fff(a, b), t);
        BB.push_back(Eigen::Triplet<T>(v.nodes.get_index(i, n),
                                       q.nodes.get_index(j, n), m));
      }
    }
  }
  return BB;
}

template <typename T, int d, int e>
std::vector<Eigen::Triplet<T>>
build_edges(const femib::finite_element_space::finite_element_space<T, d, e> &s,
            const std::function<T(femib::types::dvec<T, d>)> &b) {
  std::vector<Eigen::Triplet<T>> B;
  for (int i : s.nodes.E) {
    B.push_back(Eigen::Triplet<T>(i, 0, b(s.nodes.P[i])));
  }
  return B;
}

template <typename T, int d>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> build_zero_mean_edges(
    const femib::finite_element_space::finite_element_space<T, d, 1> &v,
    const femib::gauss::rule<T, d> &rule) {
  std::vector<Eigen::Triplet<T>> B;
  for (int n = 0; n < v.mesh.T.size(); ++n) {
    femib::types::dtrian<T, d> t = v.mesh[n];
    for (int i = 0; i < v.finite_element.base_functions.size(); ++i) {
      femib::types::F<T, d, 1> a =
          femib::util::base_function2real_function<T, d, 1>(v, n, i);

      auto g = [&](const femib::types::dvec<T, d> &x) { return a.x(x)(0); };

      T m = femib::mesh::integrate<T, d>(rule, g, t);
      B.push_back(Eigen::Triplet<T>(0, v.nodes.get_index(i, n), m));
    }
  }

  Eigen::Matrix<T, 1, Eigen::Dynamic> BB =
      femib::util::triplets2dense<T>(B, 1, v.nodes.P.size());

  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> BBB =
      Eigen::ArrayXXf::Zero(v.nodes.P.size(), v.nodes.P.size());

  T b1 = BB(0);
  BBB(0, 0) = -1;
  for (int i = 1; i < BB.cols(); ++i) {
    BBB(0, i) = BB(i) / b1;
  }
  return BBB;
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

} // namespace femib::util
#endif
