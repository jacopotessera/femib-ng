#ifndef FEMIB_NAVIER_STOKES_HPP_INCLUDED_
#define FEMIB_NAVIER_STOKES_HPP_INCLUDED_

#include "../femib/femib.hpp"
#include "../finite_element_space/finite_element_space.hpp"
#include "../gauss/gauss.hpp"
#include "navier_stokes_common.hpp"
#include "nonlinear_solvers.hpp"
#include "sparse_solvers.hpp"
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
  std::unique_ptr<femib::util::nonlinear_solver<T>> nonlinear_solver;

  navier_stokes(
      femib::finite_element_space::finite_element_space<T, d, d>
          velocity_space_,
      femib::finite_element_space::finite_element_space<T, d, 1>
          pressure_space_,
      femib::gauss::rule<T, d> quadrature_, T density_, T viscosity_,
      T time_step_,
      std::function<femib::types::dvec<T, d>(femib::types::dvec<T, d>, T)>
          body_force_ = femib::stokes_steady::default_force<T, d>(),
      std::unique_ptr<femib::util::nonlinear_solver<T>> nonlinear_solver_ =
          femib::util::default_nonlinear_solver_factory<T>(),
      std::unique_ptr<femib::util::sparse_solver<T>> sparse_solver_ =
          femib::util::default_sparse_solver_factory<T>())
      : femib::stokes::stokes<T, d>(
            std::move(velocity_space_), std::move(pressure_space_),
            std::move(quadrature_), density_, viscosity_, time_step_,
            std::move(body_force_), std::move(sparse_solver_)),
        nonlinear_solver(std::move(nonlinear_solver_)) {}
};

// backward-Euler timestep of the time-dependent Navier-Stokes system,
// with the convection term handled by s.nonlinear_solver
template <typename T, int d>
void advance(navier_stokes<T, d> &s,
             std::optional<Eigen::Matrix<T, Eigen::Dynamic, 1>>
                 extra_velocity_rhs = std::nullopt) {
  using vector_t = Eigen::Matrix<T, Eigen::Dynamic, 1>;
  s.time += s.time_step;

  vector_t previous_step_velocity;
  if (s.solutions.size() == 0)
    previous_step_velocity = vector_t::Zero(s.velocity_space.size(), 1);
  else
    previous_step_velocity =
        s.solutions[s.solutions.size() - 1].topRows(s.velocity_space.size());
  vector_t previous_step_inertia =
      (1 / s.time_step) * (s.mass_matrix * previous_step_velocity);
  Eigen::SparseMatrix<T> velocity_block_without_convection =
      s.stiffness_block + (1 / s.time_step) * s.mass_matrix;

  // the load vector can't change within a timestep
  femib::stokes::rebuild_velocity_part_rhs<T, d>(s, previous_step_inertia,
                                                 extra_velocity_rhs);

  auto step = [&s, &velocity_block_without_convection](
                  const vector_t &picard_velocity_dofs) -> vector_t {
    s.saddle_point_matrix =
        femib::stokes_steady::assemble_saddle_point_matrix<T>(
            velocity_block_without_convection +
                femib::navier_stokes_common::assemble_convection<T, d>(
                    s.velocity_space, s.quadrature, picard_velocity_dofs,
                    s.density),
            s.divergence_block, s.velocity_space.size(),
            s.pressure_space.size());
    femib::stokes_steady::rebuild_system<T, d>(s);
    return femib::stokes_steady::solve<T, d, 1>(s);
  };
  auto velocity = [&s](const vector_t &full) -> vector_t {
    return full.topRows(s.velocity_space.size());
  };

  s.solutions.emplace_back(
      s.nonlinear_solver->solve(previous_step_velocity, step, velocity));
}

} // namespace femib::navier_stokes
#endif
