#ifndef FEMIB_TYPES_HPP
#define FEMIB_TYPES_HPP

#include <Eigen/Core>
#include <vector>

namespace femib::types {

template <typename T, int d> using dvec = Eigen::Matrix<T, d, 1>;
template <typename T, int d> using dmat = Eigen::Matrix<T, d, d>;
template <typename T, int d> using dtrian = std::vector<dvec<T, d>>;

} // namespace femib::types
#endif
