#ifndef FEMIB_IB_HPP_INCLUDED_
#define FEMIB_IB_HPP_INCLUDED_

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

  Eigen::Matrix<T, Eigen::Dynamic, 1> extra_rhs =
      spread_force<T, d>(p.fluid.V, X, F, p.structure.dS);

  femib::stokes_t::advance<T, d>(p.fluid, extra_rhs);

  std::vector<femib::types::dvec<T, d>> U =
      interpolate_velocity<T, d>(p.fluid.V, p.fluid.solution.back(), X);

  for (size_t k = 0; k < p.structure.X.size(); ++k) {
    p.structure.X[k] += p.fluid.deltat * U[k];
  }
}

} // namespace femib::ib
#endif
