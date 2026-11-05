#ifndef FEMIB_IB_HPP_INCLUDED_
#define FEMIB_IB_HPP_INCLUDED_

#include "../femib/navier_stokes.hpp"
#include "../femib/nonlinear_solvers.hpp"
#include "../femib/sparse_solvers.hpp"
#include "../femib/stokes.hpp"
#include "coupling.hpp"
#include "structure.hpp"
#include <functional>
#include <memory>
#include <ranges>
#include <stdexcept>
#include <utility>

namespace femib::ib {

template <typename T, int d> struct ib_problem {
  femib::navier_stokes::navier_stokes<T, d>
      fluid; // TODO this can be stokes or navier_stokes, both have advance
  ring<T, d> structure;

  ib_problem(
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
      : fluid(std::move(velocity_space_), std::move(pressure_space_),
              std::move(quadrature_), density_, viscosity_, time_step_,
              std::move(body_force_), std::move(nonlinear_solver_),
              std::move(sparse_solver_)) {}
};

template <typename T, int d>
void advance_common(
    ib_problem<T, d> &p,
    const std::function<void(const Eigen::Matrix<T, Eigen::Dynamic, 1> &)>
        &fluid_step) {
  if (!(p.structure.reference_arc_length > 0)) {
    throw std::invalid_argument(
        "femib::ib::advance_common: p.structure.reference_arc_length must be "
        "positive");
  }

  std::vector<femib::types::dvec<T, d>> previous_step_positions =
      p.structure.positions;
  std::vector<femib::types::dvec<T, d>> force_densities =
      elastic_force<T, d>(p.structure) |
      std::views::transform([&p](const femib::types::dvec<T, d> &force) {
        return femib::types::dvec<T, d>(force /
                                        p.structure.reference_arc_length);
      }) |
      std::ranges::to<std::vector>();
  Eigen::Matrix<T, Eigen::Dynamic, 1> extra_rhs =
      spread_force<T, d>(p.fluid.velocity_space, previous_step_positions,
                         force_densities, p.structure.reference_arc_length);

  fluid_step(extra_rhs);

  std::vector<femib::types::dvec<T, d>> velocities = interpolate_velocity<T, d>(
      p.fluid.velocity_space, p.fluid.solutions.back(),
      previous_step_positions);

  advance_ring(p.structure, p.fluid.time_step, velocities);
}

// TODO eh
template <typename T, int d> void advance(ib_problem<T, d> &p) {
  advance_common<T, d>(
      p, [&p](const Eigen::Matrix<T, Eigen::Dynamic, 1> &extra_rhs) {
        femib::stokes::advance<T, d>(p.fluid, extra_rhs);
      });
}

// TODO eh
template <typename T, int d> void advance_navier_stokes(ib_problem<T, d> &p) {
  advance_common<T, d>(
      p, [&p](const Eigen::Matrix<T, Eigen::Dynamic, 1> &extra_rhs) {
        femib::navier_stokes::advance<T, d>(p.fluid, extra_rhs);
      });
}

} // namespace femib::ib
#endif
