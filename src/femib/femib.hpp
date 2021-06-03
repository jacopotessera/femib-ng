#ifndef FEMIB_HPP_INCLUDED_
#define FEMIB_HPP_INCLUDED_

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

namespace femib::util {

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
triplets2dense(std::vector<Eigen::Triplet<T>> M, int rows, int cols) {

  Eigen::SparseMatrix<T> sM = Eigen::SparseMatrix<T>(rows, cols);
  sM.setFromTriplets(M.begin(), M.end());
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> dM =
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(sM);
  return dM;
}

} // namespace femib::util
#endif
