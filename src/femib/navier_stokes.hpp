#ifndef FEMIB_NAVIER_STOKES_HPP_INCLUDED_
#define FEMIB_NAVIER_STOKES_HPP_INCLUDED_

#include "../gauss/gauss.hpp"
#include "navier_stokes_common.hpp"
#include "stokes.hpp"
#include "stokes_steady.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <execution>
#include <optional>
#include <stdexcept>

namespace femib::navier_stokes {

// backward-Euler timestep of the time-dependent Navier-Stokes system,
// with an inner Picard loop
template <typename T, int d>
void advance(femib::stokes::stokes<T, d> &s,
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
        A_base + femib::navier_stokes_common::assemble_convection<T, d>(
                     s.V, rule, xx, s.rho);
    s.AA = femib::stokes_steady::assemble_saddle_point_matrix<T>(
        top_left, s.B, s.V.nodes.P.size(), s.Q.nodes.P.size());

    femib::stokes::rebuild_rhs<T, d>(s, dd, extra_velocity_rhs);
    femib::stokes_steady::rebuild_system<T, d>(s);
    xx_new_full = femib::stokes_steady::solve<T, d, 1>(s);

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
