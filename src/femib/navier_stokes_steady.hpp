#ifndef FEMIB_NAVIER_STOKES_STEADY_HPP_INCLUDED_
#define FEMIB_NAVIER_STOKES_STEADY_HPP_INCLUDED_

#include "../finite_element_space/finite_element_space.hpp"
#include "../gauss/gauss.hpp"
#include "navier_stokes_common.hpp"
#include "nonlinear_solvers.hpp"
#include "sparse_solvers.hpp"
#include "stokes_steady.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <functional>
#include <memory>
#include <utility>

namespace femib::navier_stokes_steady {

template <typename T, int d>
struct navier_stokes : public femib::stokes_steady::stokes<T, d> {
  std::unique_ptr<femib::util::nonlinear_solver<T>> solver;

  navier_stokes(
      femib::finite_element_space::finite_element_space<T, d, d> v,
      femib::finite_element_space::finite_element_space<T, d, 1> q,
      const femib::gauss::rule<T, d> &rule, T rho = 1.0, T mu = 1.0,
      std::function<femib::types::dvec<T, d>(femib::types::dvec<T, d>, T)>
          force = femib::stokes_steady::default_force<T, d>(),
      std::unique_ptr<femib::util::nonlinear_solver<T>> solver =
          femib::util::default_nonlinear_solver_factory<T>(),
      std::unique_ptr<femib::util::sparse_solver<T>> linear_solver =
          femib::util::default_solver_factory<T>())
      : femib::stokes_steady::stokes<T, d>(std::move(v), std::move(q), rule,
                                           rho, mu, std::move(force),
                                           std::move(linear_solver)),
        solver(std::move(solver)) {}
};

template <typename T, int d>
Eigen::Matrix<T, Eigen::Dynamic, 1>
solve(navier_stokes<T, d> &s, const femib::gauss::rule<T, d> &rule) {
  using vector_t = Eigen::Matrix<T, Eigen::Dynamic, 1>;
  vector_t stokes_solution = femib::stokes_steady::solve<T, d, 1>(s);

  auto step = [&s, &rule](const vector_t &picard_velocity_dofs) -> vector_t {
    Eigen::SparseMatrix<T> top_left =
        s.stiffness_block +
        femib::navier_stokes_common::assemble_convection<T, d>(
            s.velocity_space, rule, picard_velocity_dofs, s.rho);
    s.saddle_point_matrix =
        femib::stokes_steady::assemble_saddle_point_matrix<T>(
            top_left, s.divergence_block, s.velocity_space.size(),
            s.pressure_space.size());
    femib::stokes_steady::rebuild_system<T, d>(s);
    return femib::stokes_steady::solve<T, d, 1>(s);
  };
  auto velocity = [&s](const vector_t &full) -> vector_t {
    return full.topRows(s.velocity_space.size());
  };

  vector_t xx = s.solver->solve(velocity(stokes_solution), step, velocity);
  s.solutions = {xx};
  return xx;
}

} // namespace femib::navier_stokes_steady
#endif
