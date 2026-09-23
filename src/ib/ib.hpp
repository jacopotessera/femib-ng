#ifndef FEMIB_IB_HPP_INCLUDED_
#define FEMIB_IB_HPP_INCLUDED_

#include "../femib/navier_stokes.hpp"
#include "../femib/stokes.hpp"
#include "coupling.hpp"
#include "structure.hpp"
#include <functional>
#include <stdexcept>
#include <utility>

namespace femib::ib {

template <typename T, int d> struct ib_problem {
  femib::stokes::stokes<T, d> fluid;
  ring<T, d> structure;

  ib_problem(
      femib::finite_element_space::finite_element_space<T, d, d> v,
      femib::finite_element_space::finite_element_space<T, d, 1> q,
      femib::gauss::rule<T, d> rule, T rho = 1.0, T mu = 1.0, T deltat = 0.1,
      std::function<femib::types::dvec<T, d>(femib::types::dvec<T, d>, T)>
          force = femib::stokes_steady::default_force<T, d>())
      : fluid(std::move(v), std::move(q), std::move(rule), rho, mu, deltat,
              std::move(force)) {}
};

template <typename T, int d>
void advance_common(
    ib_problem<T, d> &p,
    const std::function<void(const Eigen::Matrix<T, Eigen::Dynamic, 1> &)>
        &fluid_step) {
  if (!(p.structure.dS > 0)) {
    throw std::invalid_argument(
        "femib::ib::advance_common: p.structure.dS must be positive");
  }

  std::vector<femib::types::dvec<T, d>> X =
      p.structure.X; // pre-advance positions
  std::vector<femib::types::dvec<T, d>> F = elastic_force<T, d>(p.structure);

  // spread_force expects a force density
  for (auto &f : F) {
    f /= p.structure.dS;
  }
  Eigen::Matrix<T, Eigen::Dynamic, 1> extra_rhs =
      spread_force<T, d>(p.fluid.V, X, F, p.structure.dS);

  fluid_step(extra_rhs);

  std::vector<femib::types::dvec<T, d>> U =
      interpolate_velocity<T, d>(p.fluid.V, p.fluid.solution.back(), X);

  advance_ring(p.structure, p.fluid.deltat, U);
}

template <typename T, int d> void advance(ib_problem<T, d> &p) {
  advance_common<T, d>(
      p, [&p](const Eigen::Matrix<T, Eigen::Dynamic, 1> &extra_rhs) {
        femib::stokes::advance<T, d>(p.fluid, extra_rhs);
      });
}

template <typename T, int d>
void advance_navier_stokes(ib_problem<T, d> &p,
                           const femib::gauss::rule<T, d> &rule,
                           int max_picard_iters, T tol) {
  advance_common<T, d>(
      p, [&p, &rule, max_picard_iters,
          tol](const Eigen::Matrix<T, Eigen::Dynamic, 1> &extra_rhs) {
        femib::navier_stokes::advance<T, d>(p.fluid, rule, max_picard_iters,
                                            tol, extra_rhs);
      });
}

} // namespace femib::ib
#endif
