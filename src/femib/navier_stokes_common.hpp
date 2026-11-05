#ifndef FEMIB_NAVIER_STOKES_COMMON_HPP_INCLUDED_
#define FEMIB_NAVIER_STOKES_COMMON_HPP_INCLUDED_

#include "../affine/affine.hpp"
#include "../femib/femib.hpp"
#include "../finite_element_space/finite_element_space.hpp"
#include "../gauss/gauss.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <execution>
#include <numeric>
#include <vector>

namespace femib::navier_stokes_common {

// Assembles the convection matrix linearized around a given velocity field.
template <typename T, int d>
Eigen::SparseMatrix<T> assemble_convection(
    const femib::finite_element_space::finite_element_space<T, d, d>
        &velocity_space,
    const femib::gauss::rule<T, d> &rule,
    const Eigen::Matrix<T, Eigen::Dynamic, 1> &picard_velocity_dofs,
    T density) {
  const int basis_count =
      static_cast<int>(velocity_space.finite_element.base_functions.size());
  const int n_q = static_cast<int>(rule.nodes.size());
  const int n_tri = static_cast<int>(velocity_space.mesh.T.size());

  std::vector<Eigen::Triplet<T>> triplets(static_cast<size_t>(n_tri) *
                                          basis_count * basis_count);

  std::vector<int> tri_indices(n_tri);
  std::iota(tri_indices.begin(), tri_indices.end(), 0);

  std::for_each(
      std::execution::par_unseq, tri_indices.begin(), tri_indices.end(),
      [&](int n) {
        const femib::types::dtrian<T, d> &t = velocity_space.mesh[n];
        femib::types::dmat<T, d> inverse_jacobian =
            femib::affine::inverse_jacobian(t);
        T abs_jacobian_determinant = femib::affine::abs_jacobian_determinant(t);

        std::vector<femib::types::dvec<T, d>> phi(
            static_cast<size_t>(basis_count) * n_q);
        std::vector<femib::types::rmat<T, d, d>> dphi(
            static_cast<size_t>(basis_count) * n_q);
        std::vector<femib::types::dvec<T, d>> w_at_q(
            n_q, femib::types::dvec<T, d>::Zero());

        for (int q = 0; q < n_q; ++q) {
          const femib::types::dvec<T, d> &x_ref =
              rule.nodes[q].quadrature_point;
          for (int k = 0; k < basis_count; ++k) {
            size_t idx = static_cast<size_t>(k) * n_q + q;
            phi[idx] =
                velocity_space.finite_element.base_functions[k].value(x_ref);
            dphi[idx] =
                inverse_jacobian.transpose() *
                velocity_space.finite_element.base_functions[k].gradient(x_ref);
            int global_k = velocity_space.nodes.global_index(k, n);
            w_at_q[q] += picard_velocity_dofs(global_k) * phi[idx];
          }
        }

        for (int i = 0; i < basis_count; ++i) {
          for (int j = 0; j < basis_count; ++j) {
            T local_entry = T(0);
            for (int q = 0; q < n_q; ++q) {
              femib::types::dvec<T, d> conv_j =
                  dphi[static_cast<size_t>(j) * n_q + q].transpose() *
                  w_at_q[q];
              local_entry += rule.nodes[q].weight * abs_jacobian_determinant *
                             phi[static_cast<size_t>(i) * n_q + q].dot(conv_j);
            }
            size_t idx =
                (static_cast<size_t>(n) * basis_count + i) * basis_count + j;
            triplets[idx] = Eigen::Triplet<T>(
                velocity_space.nodes.global_index(i, n),
                velocity_space.nodes.global_index(j, n), local_entry);
          }
        }
      });

  return density * femib::util::triplets2sparse(triplets,
                                                velocity_space.nodes.P.size(),
                                                velocity_space.nodes.P.size());
}

} // namespace femib::navier_stokes_common
#endif
