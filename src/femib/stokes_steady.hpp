#ifndef FEMIB_STOKES_STEADY_HPP_INCLUDED_
#define FEMIB_STOKES_STEADY_HPP_INCLUDED_

#include "../femib/femib.hpp"
#include "../finite_element_space/finite_element_space.hpp"
#include "../gauss/gauss.hpp"
#include "../types/operators.hpp"
#include "sparse_solvers.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <functional>
#include <memory>
#include <spdlog/spdlog.h>
#include <utility>
#include <vector>

namespace femib::stokes_steady {

template <typename T, int d>
std::function<femib::types::dvec<T, d>(femib::types::dvec<T, d>, T)>
default_force() {
  return [](femib::types::dvec<T, d>, T) -> femib::types::dvec<T, d> {
    return femib::types::dvec<T, d>::Zero();
  };
}

template <typename T, int d>
std::function<T(femib::types::dvec<T, d>)>
viscous_form(femib::types::function_with_gradient<T, d, d> u,
             femib::types::function_with_gradient<T, d, d> v, T mu) {
  return [u, v, mu](const femib::types::dvec<T, d> &x) {
    return T(2.0) * mu *
           femib::operators::symmetric_gradient_double_dot(u, v)(x);
  };
}

template <typename T, int d>
std::function<T(femib::types::dvec<T, d>)>
divergence_form(femib::types::function_with_gradient<T, d, d> u,
                femib::types::function_with_gradient<T, d, 1> q) {
  return [u, q](const femib::types::dvec<T, d> &x) {
    return T(-1.0) * femib::operators::divergence(u)(x) * q.value(x)(0);
  };
}

// force at time t tested against the base function a
template <typename T, int d>
std::function<T(femib::types::dvec<T, d>)> external_force(
    femib::types::function_with_gradient<T, d, d> a,
    const std::function<femib::types::dvec<T, d>(femib::types::dvec<T, d>)>
        &force_at_t) {
  return [a, force_at_t](const femib::types::dvec<T, d> &x) {
    femib::types::dvec<T, d> fx = force_at_t(x);
    T sum = 0;
    for (int k = 0; k < d; ++k) {
      sum += fx(k) * a.value(x)(k); // TODO ?
    }
    return sum;
  };
}

template <typename T>
Eigen::SparseMatrix<T>
assemble_saddle_point_matrix(const Eigen::SparseMatrix<T> &stiffness_block,
                             const Eigen::SparseMatrix<T> &divergence_block,
                             const int nV, const int nQ) {
  std::vector<Eigen::Triplet<T>> triplets;
  triplets.reserve(static_cast<size_t>(stiffness_block.nonZeros()) +
                   2 * static_cast<size_t>(divergence_block.nonZeros()));
  for (int k = 0; k < stiffness_block.outerSize(); ++k) {
    for (typename Eigen::SparseMatrix<T>::InnerIterator it(stiffness_block, k);
         it; ++it) {
      triplets.push_back(Eigen::Triplet<T>(it.row(), it.col(), it.value()));
    }
  }
  for (int k = 0; k < divergence_block.outerSize(); ++k) {
    for (typename Eigen::SparseMatrix<T>::InnerIterator it(divergence_block, k);
         it; ++it) {
      triplets.push_back(Eigen::Triplet<T>(static_cast<int>(it.row()),
                                           nV + static_cast<int>(it.col()),
                                           it.value()));
      triplets.push_back(Eigen::Triplet<T>(nV + static_cast<int>(it.col()),
                                           static_cast<int>(it.row()),
                                           it.value()));
    }
  }
  Eigen::SparseMatrix<T> saddle_point_matrix(nV + nQ, nV + nQ);
  saddle_point_matrix.setFromTriplets(triplets.begin(), triplets.end());
  return saddle_point_matrix;
}

// entry(j) = integral of pressure basis function j over the whole domain, so
// entry . q = integral of q over the domain: the zero-mean pressure gauge
// constraint. Equivalent to build_load_vector(v, rule, load≡1). Un-normalized
// (not divided by domain volume), which is fine since the constraint is
// homogeneous (scaling it by a positive constant doesn't change the solution).
template <typename T, int d>
Eigen::Matrix<T, Eigen::Dynamic, 1> build_pressure_constraint_vector(
    const femib::finite_element_space::finite_element_space<T, d, 1>
        &v, // TODO no need to pass finite_element_space, it's the pressure one
    const femib::gauss::rule<T, d> &rule) {
  return femib::util::build_load_vector<T, d, 1>(
      v, rule, [](femib::types::function_with_gradient<T, d, 1> a) {
        return [a](const femib::types::dvec<T, d> &x) { return a.value(x)(0); };
      });
}

template <typename T, int d> struct stokes {
  femib::finite_element_space::finite_element_space<T, d, d> velocity_space;
  femib::finite_element_space::finite_element_space<T, d, 1> pressure_space;
  T density;
  T viscosity;
  std::function<femib::types::dvec<T, d>(femib::types::dvec<T, d>, T)>
      body_force;

  Eigen::SparseMatrix<T> stiffness_block;
  Eigen::SparseMatrix<T> divergence_block;
  Eigen::Matrix<T, Eigen::Dynamic, 1> boundary_values_on_nodes;
  Eigen::Matrix<T, Eigen::Dynamic, 1> pressure_constraint_vector;
  std::vector<int> not_edge_nodes;

  Eigen::SparseMatrix<T> saddle_point_matrix;
  Eigen::Matrix<T, Eigen::Dynamic, 1> rhs_vector;
  femib::util::sparse_linear_system<T> reduced_system;
  std::unique_ptr<femib::util::sparse_solver<T>> sparse_solver;

  std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1>> solutions;

  stokes(femib::finite_element_space::finite_element_space<T, d, d>
             velocity_space_,
         femib::finite_element_space::finite_element_space<T, d, 1>
             pressure_space_,
         const femib::gauss::rule<T, d> &integration_rule, T density_,
         T viscosity_,
         std::function<femib::types::dvec<T, d>(femib::types::dvec<T, d>, T)>
             body_force_ = default_force<T, d>(),
         std::unique_ptr<femib::util::sparse_solver<T>> sparse_solver_ =
             femib::util::default_sparse_solver_factory<T>());

  std::vector<std::pair<types::dvec<T, d>, types::dvec<T, d>>>
  plot_velocity(T delta) {
    return velocity_space.plot(solutions.back().topRows(velocity_space.size()),
                               delta);
  }
  std::vector<std::pair<types::dvec<T, d>, types::dvec<T, 1>>>
  plot_pressure(T delta) {
    return pressure_space.plot(solutions.back().tail(pressure_space.size()),
                               delta);
  }
};

template <typename T, int d>
femib::util::sparse_linear_system<T>
augment_with_pressure_gauge(const stokes<T, d> &s) {

  femib::util::sparse_linear_system<T> base =
      femib::util::build_reduced_system<T>(s.saddle_point_matrix, s.rhs_vector,
                                           s.boundary_values_on_nodes,
                                           s.not_edge_nodes);

  int n = static_cast<int>(s.not_edge_nodes.size());
  Eigen::Matrix<T, Eigen::Dynamic, 1> constraint_row_reduced =
      Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(n);
  for (int k = 0; k < n; ++k) {
    int global_i = s.not_edge_nodes[k];
    if (global_i >= s.velocity_space.size()) {
      int pressure_i = global_i - s.velocity_space.size();
      constraint_row_reduced(k) = s.pressure_constraint_vector(pressure_i);
    }
  }

  std::vector<Eigen::Triplet<T>> triplets;
  triplets.reserve(base.A.nonZeros() + 2 * n + (n + 1));
  for (int k = 0; k < base.A.outerSize(); ++k) {
    for (typename Eigen::SparseMatrix<T>::InnerIterator it(base.A, k); it;
         ++it) {
      triplets.push_back(Eigen::Triplet<T>(it.row(), it.col(), it.value()));
    }
  }
  for (int k = 0; k < n; ++k) {
    if (constraint_row_reduced(k) != T(0)) {
      triplets.push_back(Eigen::Triplet<T>(k, n, constraint_row_reduced(k)));
      triplets.push_back(Eigen::Triplet<T>(n, k, constraint_row_reduced(k)));
    }
  }

  // TODO this must be a "regularize problem" , tikhonov as implementation

  // Tikhonov regularization
  T reg = T(1e-8) * s.stiffness_block.diagonal().cwiseAbs().maxCoeff();
  if (reg == T(0)) {
    SPDLOG_WARN("femib::stokes_steady::augment_with_pressure_gauge: "
                "s.stiffness_block's "
                "diagonal is all-zero: system may be left singular");
  }
  int nV_total = s.velocity_space.size();
  for (int k = 0; k < n; ++k) {
    bool is_velocity = s.not_edge_nodes[k] < nV_total;
    triplets.push_back(Eigen::Triplet<T>(k, k, is_velocity ? reg : -reg));
  }
  triplets.push_back(Eigen::Triplet<T>(n, n, -reg));

  Eigen::SparseMatrix<T> AAA_aug(n + 1, n + 1);
  AAA_aug.setFromTriplets(triplets.begin(), triplets.end());

  Eigen::Matrix<T, Eigen::Dynamic, 1> bbb_aug =
      Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(n + 1);
  bbb_aug.topRows(n) = base.b;

  return {AAA_aug, bbb_aug};
}

// TODO for not steady problems?
template <typename T, int d> void rebuild_system(stokes<T, d> &s) {
  s.reduced_system = augment_with_pressure_gauge<T, d>(s);
}

template <typename T, int d, int e>
Eigen::Matrix<T, Eigen::Dynamic, 1> solve(const stokes<T, d> &s) {

  Eigen::Matrix<T, Eigen::Dynamic, 1> x =
      s.sparse_solver->solve(s.reduced_system.A, s.reduced_system.b);

  // drop the last row, constraint on pressure
  Eigen::Matrix<T, Eigen::Dynamic, 1> x_no_lambda = x.topRows(x.rows() - 1);

  return femib::util::build_full_solution<T>(
      x_no_lambda, s.boundary_values_on_nodes,
      s.velocity_space.size() + s.pressure_space.size(), s.not_edge_nodes,
      s.velocity_space.nodes.E);
}

template <typename T, int d>
stokes<T, d>::stokes(
    femib::finite_element_space::finite_element_space<T, d, d> velocity_space_,
    femib::finite_element_space::finite_element_space<T, d, 1> pressure_space_,
    const femib::gauss::rule<T, d> &integration_rule, T density_, T viscosity_,
    std::function<femib::types::dvec<T, d>(femib::types::dvec<T, d>, T)>
        body_force_,
    std::unique_ptr<femib::util::sparse_solver<T>> sparse_solver_)
    : velocity_space(std::move(velocity_space_)),
      pressure_space(std::move(pressure_space_)), density(density_),
      viscosity(viscosity_), body_force(std::move(body_force_)),
      sparse_solver(std::move(sparse_solver_)) {
  stiffness_block = femib::util::triplets2sparse(
      femib::util::build_diagonal_block<T, d, d>(
          velocity_space, integration_rule,
          [this](femib::types::function_with_gradient<T, d, d> u,
                 femib::types::function_with_gradient<T, d, d> w) {
            return viscous_form<T, d>(u, w, viscosity);
          }),
      velocity_space.size(), velocity_space.size());
  divergence_block = femib::util::triplets2sparse(
      femib::util::build_off_diagonal_block<T, d>(
          velocity_space, pressure_space, integration_rule,
          divergence_form<T, d>),
      velocity_space.size(), pressure_space.size());

  std::function<femib::types::dvec<T, d>(femib::types::dvec<T, d>)> force_at_0 =
      [f = this->body_force](const femib::types::dvec<T, d> &x) {
        return f(x, T(0));
      };
  rhs_vector = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(
      velocity_space.size() + pressure_space.size(), 1);
  rhs_vector.block(0, 0, velocity_space.size(), 1) =
      femib::util::build_load_vector<T, d, d>(
          velocity_space, integration_rule,
          [force_at_0](femib::types::function_with_gradient<T, d, d> a) {
            return external_force<T, d>(a, force_at_0);
          });

  std::function<T(femib::types::dvec<T, d>)> zero_boundary =
      [](const femib::types::dvec<T, d> &) { return T(0); };
  boundary_values_on_nodes = femib::util::triplets2dense(
      femib::util::build_boundary_values_on_nodes<T, d, d>(velocity_space,
                                                           zero_boundary),
      velocity_space.size() + pressure_space.size(), 1);
  pressure_constraint_vector =
      build_pressure_constraint_vector<T, d>(pressure_space, integration_rule);
  not_edge_nodes = femib::util::find_not_edge_nodes<T, d, d>(velocity_space);
  for (int i = 0; i < pressure_space.size(); ++i) {
    not_edge_nodes.push_back(velocity_space.size() + i);
  }
  saddle_point_matrix = assemble_saddle_point_matrix<T>(
      stiffness_block, divergence_block, velocity_space.size(),
      pressure_space.size());
  rebuild_system<T, d>(*this);
}

} // namespace femib::stokes_steady
#endif
