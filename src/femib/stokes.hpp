#ifndef FEMIB_STOKES_HPP_INCLUDED_
#define FEMIB_STOKES_HPP_INCLUDED_

#include "../femib/femib.hpp"
#include "../finite_element_space/finite_element_space.hpp"
#include "../gauss/gauss.hpp"
#include "stokes_steady.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <functional>
#include <optional>
#include <utility>

namespace femib::stokes {

template <typename T, int d>
struct stokes : public femib::stokes_steady::stokes<T, d> {
  femib::gauss::rule<T, d> rule;
  T deltat;
  Eigen::SparseMatrix<T> M;
  T time = 0;

  stokes(femib::finite_element_space::finite_element_space<T, d, d> v,
         femib::finite_element_space::finite_element_space<T, d, 1> q,
         femib::gauss::rule<T, d> rule, T rho = 1.0, T mu = 1.0,
         T deltat = 0.1, // TODO eh
         std::function<femib::types::dvec<T, d>(femib::types::dvec<T, d>, T)>
             force = femib::stokes_steady::default_force<T, d>())
      : femib::stokes_steady::stokes<T, d>(std::move(v), std::move(q), rule,
                                           rho, mu, std::move(force)),
        rule(std::move(rule)), deltat(deltat) {
    M = this->rho *
        femib::util::build_mass_matrix<T, d, d>(this->V, this->rule);
  }
};

// TODO we need to give better names to stuff...
// Sets ff's velocity block for the step ending at s.time (the caller must
// already have advanced s.time): s.force at s.time tested against V, plus dd
// (the backward-Euler mass term's known part), plus extra_velocity_rhs.
template <typename T, int d>
void rebuild_rhs(
    stokes<T, d> &s, const Eigen::Matrix<T, Eigen::Dynamic, 1> &dd,
    std::optional<Eigen::Matrix<T, Eigen::Dynamic, 1>> extra_velocity_rhs) {

  std::function<femib::types::dvec<T, d>(femib::types::dvec<T, d>)> force_n1 =
      [force = s.force, time_n1 = s.time](const femib::types::dvec<T, d> &x) {
        return force(x, time_n1);
      };
  auto ggg = [force_n1](femib::types::F<T, d, d> a) {
    return femib::stokes_steady::external_force<T, d>(a, force_n1);
  };
  Eigen::Matrix<T, Eigen::Dynamic, 1> velocity_rhs =
      femib::util::build_load_vector<T, d, d>(s.V, s.rule, ggg) + dd;

  if (extra_velocity_rhs.has_value()) {
    velocity_rhs += extra_velocity_rhs.value();
  }

  s.ff.block(0, 0, s.V.size(), 1) = velocity_rhs;
}

template <typename T, int d>
void advance(stokes<T, d> &s, std::optional<Eigen::Matrix<T, Eigen::Dynamic, 1>>
                                  extra_velocity_rhs = std::nullopt) {
  s.time += s.deltat;

  Eigen::Matrix<T, Eigen::Dynamic, 1> u_1;
  if (s.solution.size() == 0)
    u_1 = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(s.V.size(), 1);
  else
    // last timestep velocity
    u_1 = s.solution[s.solution.size() - 1].topRows(s.V.size());

  Eigen::SparseMatrix<T> DD = (1 / s.deltat) * s.M;
  Eigen::Matrix<T, Eigen::Dynamic, 1> dd = (1 / s.deltat) * (s.M * u_1);

  int nV = s.V.size();
  int nQ = s.Q.size();
  Eigen::SparseMatrix<T> top_left = s.A + DD;
  s.AA = femib::stokes_steady::assemble_saddle_point_matrix<T>(top_left, s.B,
                                                               nV, nQ);

  femib::stokes::rebuild_rhs<T, d>(s, dd, extra_velocity_rhs);
  femib::stokes_steady::rebuild_system<T, d>(s);

  Eigen::Matrix<T, Eigen::Dynamic, 1> xx =
      femib::stokes_steady::solve<T, d, 1>(s);

  s.solution.emplace_back(xx);
}

} // namespace femib::stokes
#endif
