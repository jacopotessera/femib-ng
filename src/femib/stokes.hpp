#ifndef FEMIB_STOKES_HPP_INCLUDED_
#define FEMIB_STOKES_HPP_INCLUDED_

#include "../femib/femib.hpp"
#include "../finite_element_space/finite_element_space.hpp"
#include "../gauss/gauss.hpp"
#include "sparse_solvers.hpp"
#include "stokes_steady.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <functional>
#include <memory>
#include <optional>
#include <utility>

namespace femib::stokes {

template <typename T, int d>
struct stokes : public femib::stokes_steady::stokes<T, d> {
  femib::gauss::rule<T, d> quadrature;
  T time_step;
  Eigen::SparseMatrix<T> M;
  T time = 0;

  stokes(femib::finite_element_space::finite_element_space<T, d, d>
             velocity_space_,
         femib::finite_element_space::finite_element_space<T, d, 1>
             pressure_space_,
         femib::gauss::rule<T, d> quadrature_, T density_ = 1.0,
         T viscosity_ = 1.0,
         T time_step_ = 0.1, // TODO eh
         std::function<femib::types::dvec<T, d>(femib::types::dvec<T, d>, T)>
             body_force_ = femib::stokes_steady::default_force<T, d>(),
         std::unique_ptr<femib::util::sparse_solver<T>> sparse_solver_ =
             femib::util::default_solver_factory<T>())
      : femib::stokes_steady::stokes<T, d>(
            std::move(velocity_space_), std::move(pressure_space_), quadrature_,
            density_, viscosity_, std::move(body_force_),
            std::move(sparse_solver_)),
        quadrature(std::move(quadrature_)), time_step(time_step_) {
    M = this->density * femib::util::build_mass_matrix<T, d, d>(
                            this->velocity_space, this->quadrature);
  }
};

// TODO we need to give better names to stuff...
// Sets ff's velocity block for the step ending at s.time (the caller must
// already have advanced s.time): s.body_force at s.time tested against V, plus
// dd (the backward-Euler mass term's known part), plus extra_velocity_rhs.
template <typename T, int d>
void rebuild_rhs(
    stokes<T, d> &s, const Eigen::Matrix<T, Eigen::Dynamic, 1> &dd,
    std::optional<Eigen::Matrix<T, Eigen::Dynamic, 1>> extra_velocity_rhs) {

  std::function<femib::types::dvec<T, d>(femib::types::dvec<T, d>)> force_n1 =
      [force = s.body_force, time_n1 = s.time](
          const femib::types::dvec<T, d> &x) { return force(x, time_n1); };
  auto ggg = [force_n1](femib::types::F<T, d, d> a) {
    return femib::stokes_steady::external_force<T, d>(a, force_n1);
  };
  Eigen::Matrix<T, Eigen::Dynamic, 1> velocity_rhs =
      femib::util::build_load_vector<T, d, d>(s.velocity_space, s.quadrature,
                                              ggg) +
      dd;

  if (extra_velocity_rhs.has_value()) {
    velocity_rhs += extra_velocity_rhs.value();
  }

  s.rhs_vector.block(0, 0, s.velocity_space.size(), 1) = velocity_rhs;
}

template <typename T, int d>
void advance(stokes<T, d> &s, std::optional<Eigen::Matrix<T, Eigen::Dynamic, 1>>
                                  extra_velocity_rhs = std::nullopt) {
  s.time += s.time_step;

  Eigen::Matrix<T, Eigen::Dynamic, 1> u_1;
  if (s.solutions.size() == 0)
    u_1 = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(s.velocity_space.size(), 1);
  else
    // last timestep velocity
    u_1 = s.solutions[s.solutions.size() - 1].topRows(s.velocity_space.size());

  Eigen::SparseMatrix<T> DD = (1 / s.time_step) * s.M;
  Eigen::Matrix<T, Eigen::Dynamic, 1> dd = (1 / s.time_step) * (s.M * u_1);

  const int nV = s.velocity_space.size();
  const int nQ = s.pressure_space.size();
  Eigen::SparseMatrix<T> top_left = s.stiffness_block + DD;
  s.saddle_point_matrix = femib::stokes_steady::assemble_saddle_point_matrix<T>(
      top_left, s.divergence_block, nV, nQ);

  femib::stokes::rebuild_rhs<T, d>(s, dd, extra_velocity_rhs);
  femib::stokes_steady::rebuild_system<T, d>(s);

  Eigen::Matrix<T, Eigen::Dynamic, 1> xx =
      femib::stokes_steady::solve<T, d, 1>(s);

  s.solutions.emplace_back(xx);
}

} // namespace femib::stokes
#endif
