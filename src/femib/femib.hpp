#ifndef FEMIB_HPP_INCLUDED_
#define FEMIB_HPP_INCLUDED_

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>

#include "../affine/affine.hpp"
#include "../finite_element_space/finite_element_space.hpp"
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

} // namespace femib::util
#endif
