#ifndef FEMIB_IB_HPP_INCLUDED_
#define FEMIB_IB_HPP_INCLUDED_

#include "../femib/navier_stokes.hpp"
#include "../femib/stokes_t.hpp"
#include "coupling.hpp"
#include "structure.hpp"

namespace femib::ib {

template <typename T, int d> struct ib_problem {
  femib::stokes_t::stokes<T, d> fluid;
  ring<T, d> structure;
};

template <typename T, int d>
void init(ib_problem<T, d> &p, const femib::gauss::rule<T, d> &rule) {
  femib::stokes_t::init<T, d>(p.fluid, rule);
}

template <typename T, int d> void advance(ib_problem<T, d> &p) {
  std::vector<femib::types::dvec<T, d>> X =
      p.structure.X; // pre-advance positions
  std::vector<femib::types::dvec<T, d>> F = elastic_force<T, d>(p.structure);

  // TODO spread_force expects a force DENSITY
  Eigen::Matrix<T, Eigen::Dynamic, 1> extra_rhs =
      spread_force<T, d>(p.fluid.V, X, F, p.structure.dS);

  femib::stokes_t::advance<T, d>(p.fluid, extra_rhs);

  std::vector<femib::types::dvec<T, d>> U =
      interpolate_velocity<T, d>(p.fluid.V, p.fluid.solution.back(), X);

  advance_ring(p.structure, p.fluid.deltat, U);
}

template <typename T, int d>
void advance_navier_stokes(ib_problem<T, d> &p,
                           const femib::gauss::rule<T, d> &rule, T reynolds,
                           int max_picard_iters, T tol) {
  std::vector<femib::types::dvec<T, d>> X = p.structure.X;
  std::vector<femib::types::dvec<T, d>> F = elastic_force<T, d>(p.structure);

  // TODO spread_force expects a force DENSITY
  Eigen::Matrix<T, Eigen::Dynamic, 1> extra_rhs =
      spread_force<T, d>(p.fluid.V, X, F, p.structure.dS);

  femib::navier_stokes::advance<T, d>(p.fluid, rule, reynolds, max_picard_iters,
                                      tol, extra_rhs);

  std::vector<femib::types::dvec<T, d>> U =
      interpolate_velocity<T, d>(p.fluid.V, p.fluid.solution.back(), X);

  advance_ring(p.structure, p.fluid.deltat, U);
}

} // namespace femib::ib
#endif
