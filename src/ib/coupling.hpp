#ifndef FEMIB_IB_COUPLING_HPP_INCLUDED_
#define FEMIB_IB_COUPLING_HPP_INCLUDED_

#include "../finite_element_space/finite_element_space.hpp"
#include "../types/types.hpp"
#include <vector>

namespace femib::ib {

template <typename T, int d>
std::vector<femib::types::dvec<T, d>> interpolate_velocity(
    femib::finite_element_space::finite_element_space<T, d, d> &V,
    const Eigen::Matrix<T, Eigen::Dynamic, 1> &u,
    std::vector<femib::types::dvec<T, d>> &points) {
  return V.interpolate(u, points);
}

template <typename T, int d>
Eigen::Matrix<T, Eigen::Dynamic, 1>
spread_force(femib::finite_element_space::finite_element_space<T, d, d> &V,
             std::vector<femib::types::dvec<T, d>> &points,
             const std::vector<femib::types::dvec<T, d>> &forces, T dS) {

  Eigen::Matrix<T, Eigen::Dynamic, 1> result =
      Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(V.nodes.P.size());

  std::vector<int> points_positions =
      femib::mesh::find_points<T, d>(V.mesh, points);

  for (size_t k = 0; k < points.size(); ++k) {
    if (points_positions[k] < 0) {
      continue;
    }
    femib::types::dtrian<T, d> t = V.mesh[points_positions[k]];
    femib::types::dvec<T, d> ref =
        femib::affine::affine_inv<T, d>(t, points[k]);
    for (int j = 0; j < V.finite_element.base_functions.size(); ++j) {
      femib::types::F<T, d, d> phi = V.finite_element.base_functions[j];
      // variational formulation of the force spreading
      T contribution = forces[k].dot(phi.x(ref)) * dS;
      result(V.nodes.get_index(j, points_positions[k])) += contribution;
    }
  }
  return result;
}

} // namespace femib::ib
#endif
