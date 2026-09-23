#ifndef FEMIB_NAVIER_STOKES_HPP_INCLUDED_
#define FEMIB_NAVIER_STOKES_HPP_INCLUDED_

#include "../femib/femib.hpp"
#include "../finite_element_space/finite_element_space.hpp"
#include "../gauss/gauss.hpp"
#include "navier_stokes_common.hpp"
#include "nonlinear_solvers.hpp"
#include "stokes.hpp"
#include "stokes_steady.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <functional>
#include <memory>
#include <optional>
#include <utility>

namespace femib::navier_stokes {

template <typename T, int d>
struct navier_stokes : public femib::stokes::stokes<T, d> {
  std::unique_ptr<femib::util::nonlinear_solver<T>> solver;

  navier_stokes(
      femib::finite_element_space::finite_element_space<T, d, d> v,
      femib::finite_element_space::finite_element_space<T, d, 1> q,
      femib::gauss::rule<T, d> rule, T rho = 1.0, T mu = 1.0, T deltat = 0.1,
      std::function<femib::types::dvec<T, d>(femib::types::dvec<T, d>, T)>
          force = femib::stokes_steady::default_force<T, d>(),
      std::unique_ptr<femib::util::nonlinear_solver<T>> solver =
          femib::util::default_nonlinear_solver_factory<T>())
      : femib::stokes::stokes<T, d>(std::move(v), std::move(q), std::move(rule),
                                    rho, mu, deltat, std::move(force)),
        solver(std::move(solver)) {}
};

// backward-Euler timestep of the time-dependent Navier-Stokes system,
// with the convection term handled by s.solver
template <typename T, int d>
void advance(navier_stokes<T, d> &s,
             std::optional<Eigen::Matrix<T, Eigen::Dynamic, 1>>
                 extra_velocity_rhs = std::nullopt) {
  using vector_t = Eigen::Matrix<T, Eigen::Dynamic, 1>;
  s.time += s.deltat;

  vector_t u_1;
  if (s.solution.size() == 0)
    u_1 = vector_t::Zero(s.V.size(), 1);
  else
    u_1 = s.solution[s.solution.size() - 1].topRows(s.V.size());
  Eigen::SparseMatrix<T> DD = (1 / s.deltat) * s.M;
  vector_t dd = (1 / s.deltat) * (s.M * u_1);
  Eigen::SparseMatrix<T> A_base = s.A + DD;

  // the load vector can't change within a timestep
  femib::stokes::rebuild_rhs<T, d>(s, dd, extra_velocity_rhs);

  auto step = [&s, &A_base](const vector_t &picard_velocity_dofs) -> vector_t {
    Eigen::SparseMatrix<T> top_left =
        A_base + femib::navier_stokes_common::assemble_convection<T, d>(
                     s.V, s.rule, picard_velocity_dofs, s.rho);
    s.AA = femib::stokes_steady::assemble_saddle_point_matrix<T>(
        top_left, s.B, s.V.size(), s.Q.size());
    femib::stokes_steady::rebuild_system<T, d>(s);
    return femib::stokes_steady::solve<T, d, 1>(s);
  };
  auto velocity = [&s](const vector_t &full) -> vector_t {
    return full.topRows(s.V.size());
  };

  s.solution.emplace_back(s.solver->solve(u_1, step, velocity));
}

} // namespace femib::navier_stokes
#endif
