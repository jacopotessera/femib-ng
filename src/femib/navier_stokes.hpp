#ifndef FEMIB_NAVIER_STOKES_HPP_INCLUDED_
#define FEMIB_NAVIER_STOKES_HPP_INCLUDED_

#include "../affine/affine.hpp"
#include "../femib/femib.hpp"
#include "../finite_element_space/finite_element_space.hpp"
#include "../gauss/gauss.hpp"
#include "stokes.hpp"
#include "stokes_steady.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <execution>
#include <numeric>
#include <optional>
#include <stdexcept>

namespace femib::navier_stokes {

// convection term, using Picard's iterative method
template <typename T, int d>
Eigen::SparseMatrix<T> assemble_convection(
    const femib::finite_element_space::finite_element_space<T, d, d> &V,
    const femib::gauss::rule<T, d> &rule,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>
        &w_dofs, // TODO uh? u_at_last_step
    T rho) {
  const int basis_count =
      static_cast<int>(V.finite_element.base_functions.size());
  const int n_q = static_cast<int>(rule.nodes.size());
  const int n_tri = static_cast<int>(V.mesh.T.size());

  std::vector<Eigen::Triplet<T>> BB(static_cast<size_t>(n_tri) * basis_count *
                                    basis_count);

  std::vector<int> tri_indices(n_tri);
  std::iota(tri_indices.begin(), tri_indices.end(), 0);

  std::for_each(
      std::execution::par, tri_indices.begin(), tri_indices.end(), [&](int n) {
        const femib::types::dtrian<T, d> &t = V.mesh[n];
        femib::types::dmat<T, d> Binv = femib::affine::affineBinv(t);
        femib::types::dvec<T, d> bb = femib::affine::affineb(t);
        T detB = femib::affine::affineBdet(t);

        std::vector<femib::types::dvec<T, d>> phi(
            static_cast<size_t>(basis_count) * n_q);
        std::vector<femib::types::rmat<T, d, d>> dphi(
            static_cast<size_t>(basis_count) * n_q);
        std::vector<femib::types::dvec<T, d>> w_at_q(
            n_q, femib::types::dvec<T, d>::Zero());

        for (int q = 0; q < n_q; ++q) {
          const femib::types::dvec<T, d> &x_ref = rule.nodes[q].node;
          for (int k = 0; k < basis_count; ++k) {
            size_t idx = static_cast<size_t>(k) * n_q + q;
            phi[idx] = V.finite_element.base_functions[k].x(x_ref);
            dphi[idx] =
                Binv.transpose() * V.finite_element.base_functions[k].dx(x_ref);
            int global_k = V.nodes.get_index(k, n);
            w_at_q[q] += w_dofs(global_k) * phi[idx];
          }
        }

        for (int i = 0; i < basis_count; ++i) {
          for (int j = 0; j < basis_count; ++j) {
            T m = T(0);
            for (int q = 0; q < n_q; ++q) {
              femib::types::dvec<T, d> conv_j =
                  dphi[static_cast<size_t>(j) * n_q + q].transpose() *
                  w_at_q[q];
              m += rule.nodes[q].weight * detB *
                   phi[static_cast<size_t>(i) * n_q + q].dot(conv_j);
            }
            size_t idx =
                (static_cast<size_t>(n) * basis_count + i) * basis_count + j;
            BB[idx] = Eigen::Triplet<T>(V.nodes.get_index(i, n),
                                        V.nodes.get_index(j, n), m);
          }
        }
      });

  return rho *
         femib::util::triplets2sparse(BB, V.nodes.P.size(), V.nodes.P.size());
}

// Picard iteration:
// solve Stokes once, then repeatedly add the convection contribution to the
// frozen Stokes matrix, rebuild the solvable system, and re-solve, until the
// relative change is below tol or max_picard_iters is reached.
template <typename T, int d>
Eigen::Matrix<T, Eigen::Dynamic, 1>
solve_steady(femib::stokes_steady::stokes<T, d> &s,
             const femib::gauss::rule<T, d> &rule, int max_picard_iters,
             T tol) {
  Eigen::SparseMatrix<T> A_stokes = s.A;
  Eigen::Matrix<T, Eigen::Dynamic, 1> xx =
      femib::stokes_steady::solve<T, d, 1>(s);
  for (int iter = 0; iter < max_picard_iters; ++iter) {
    Eigen::Matrix<T, Eigen::Dynamic, 1> w_dofs = xx.topRows(s.V.nodes.P.size());
    s.A = A_stokes + assemble_convection<T, d>(s.V, rule, w_dofs, s.rho);
    femib::stokes_steady::rebuild_system<T, d>(s, rule);
    Eigen::Matrix<T, Eigen::Dynamic, 1> xx_new =
        femib::stokes_steady::solve<T, d, 1>(s);
    T rel_change =
        (xx_new - xx).norm() / std::max(xx_new.norm(), static_cast<T>(1e-8));
    xx = xx_new;
    if (rel_change < tol)
      break;
  }
  // Restore s.A to the pure-Stokes stiffness
  s.A = A_stokes;
  return xx;
}

// backward-Euler timestep of the time-dependent Navier-Stokes system,
// with an inner Picard loop
template <typename T, int d>
void advance(femib::stokes_t::stokes<T, d> &s,
             const femib::gauss::rule<T, d> &rule, int max_picard_iters, T tol,
             std::optional<Eigen::Matrix<T, Eigen::Dynamic, 1>>
                 extra_velocity_rhs = std::nullopt) {
  s.time += s.deltat;

  if (max_picard_iters <= 0) {
    throw std::invalid_argument(
        "femib::navier_stokes::advance: max_picard_iters must be >= 1");
  }

  Eigen::Matrix<T, Eigen::Dynamic, 1> u_1;
  if (s.solution.size() == 0)
    u_1 = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(s.V.nodes.P.size(), 1);
  else
    u_1 = s.solution[s.solution.size() - 1].topRows(s.V.nodes.P.size());
  Eigen::SparseMatrix<T> DD = (1 / s.deltat) * s.M;
  Eigen::Matrix<T, Eigen::Dynamic, 1> dd = (1 / s.deltat) * (s.M * u_1);
  Eigen::SparseMatrix<T> A_base = s.A + DD;

  // Picard loop
  Eigen::Matrix<T, Eigen::Dynamic, 1> xx = u_1;
  Eigen::Matrix<T, Eigen::Dynamic, 1> xx_new_full;
  for (int iter = 0; iter < max_picard_iters; ++iter) {
    Eigen::SparseMatrix<T> top_left =
        A_base + assemble_convection<T, d>(s.V, rule, xx, s.rho);
    s.AA = femib::stokes_t::assemble_saddle_point_matrix<T>(
        top_left, s.B, s.V.nodes.P.size(), s.Q.nodes.P.size());

    femib::stokes_t::rebuild_system<T, d>(s, dd, extra_velocity_rhs);
    xx_new_full = femib::stokes_t::solve<T, d, 1>(s);

    Eigen::Matrix<T, Eigen::Dynamic, 1> xx_new_velocity =
        xx_new_full.topRows(s.V.nodes.P.size());
    T rel_change = (xx_new_velocity - xx).norm() /
                   std::max(xx_new_velocity.norm(), static_cast<T>(1e-8));
    xx = xx_new_velocity;
    if (iter > 0 && rel_change < tol) {
      break;
    }
  }
  s.solution.emplace_back(xx_new_full);
}

} // namespace femib::navier_stokes
#endif
