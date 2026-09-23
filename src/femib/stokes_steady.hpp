#ifndef FEMIB_STOKES_STEADY_HPP_INCLUDED_
#define FEMIB_STOKES_STEADY_HPP_INCLUDED_

#include "../femib/femib.hpp"
#include "../finite_element_space/finite_element_space.hpp"
#include "../gauss/gauss.hpp"
#include "../mesh/mesh.hpp"
#include "../types/differential_operation.hpp"
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <functional>
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
stokes_a(femib::types::F<T, d, d> u, femib::types::F<T, d, d> v, T mu) {
  return [u, v, mu](const femib::types::dvec<T, d> &x) {
    return T(2.0) * mu * dpi(u, v)(x);
  };
}

template <typename T, int d>
std::function<T(femib::types::dvec<T, d>)>
stokes_b(femib::types::F<T, d, d> u, femib::types::F<T, d, 1> q) {
  return [u, q](const femib::types::dvec<T, d> &x) {
    return T(-1.0) * div(u)(x) * q.x(x)(0);
  };
}

// force at time t tested against the base function a
template <typename T, int d>
std::function<T(femib::types::dvec<T, d>)> external_force(
    femib::types::F<T, d, d> a,
    const std::function<femib::types::dvec<T, d>(femib::types::dvec<T, d>)>
        &force_at_t) {
  return [a, force_at_t](const femib::types::dvec<T, d> &x) {
    femib::types::dvec<T, d> fx = force_at_t(x);
    T sum = 0;
    for (int k = 0; k < d; ++k) {
      sum += fx(k) * a.x(x)(k); // TODO ?
    }
    return sum;
  };
}

template <typename T>
Eigen::SparseMatrix<T>
assemble_saddle_point_matrix(const Eigen::SparseMatrix<T> &top_left,
                             const Eigen::SparseMatrix<T> &B, int nV, int nQ) {
  std::vector<Eigen::Triplet<T>> triplets;
  triplets.reserve((size_t)top_left.nonZeros() + 2 * (size_t)B.nonZeros());
  for (int k = 0; k < top_left.outerSize(); ++k) {
    for (typename Eigen::SparseMatrix<T>::InnerIterator it(top_left, k); it;
         ++it) {
      triplets.push_back(
          Eigen::Triplet<T>((int)it.row(), (int)it.col(), it.value()));
    }
  }
  for (int k = 0; k < B.outerSize(); ++k) {
    for (typename Eigen::SparseMatrix<T>::InnerIterator it(B, k); it; ++it) {
      triplets.push_back(
          Eigen::Triplet<T>((int)it.row(), nV + (int)it.col(), it.value()));
      triplets.push_back(
          Eigen::Triplet<T>(nV + (int)it.col(), (int)it.row(), it.value()));
    }
  }
  Eigen::SparseMatrix<T> AA(nV + nQ, nV + nQ);
  AA.setFromTriplets(triplets.begin(), triplets.end());
  return AA;
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
      v, rule, [](femib::types::F<T, d, 1> a) {
        return [a](const femib::types::dvec<T, d> &x) { return a.x(x)(0); };
      });
}

template <typename T, int d> struct stokes {
  femib::finite_element_space::finite_element_space<T, d, d> V;
  femib::finite_element_space::finite_element_space<T, d, 1> Q;
  T rho = 1.0;
  T mu = 1.0;
  std::function<femib::types::dvec<T, d>(femib::types::dvec<T, d>, T)> force =
      default_force<T, d>();

  Eigen::SparseMatrix<T> A;
  Eigen::SparseMatrix<T> B;
  Eigen::Matrix<T, Eigen::Dynamic, 1> bV;
  Eigen::Matrix<T, Eigen::Dynamic, 1> pressure_constraint_vector;
  std::vector<int> not_edges;

  Eigen::SparseMatrix<T> AA;
  Eigen::Matrix<T, Eigen::Dynamic, 1> ff;
  femib::util::sparse_solvable_equations<T> solvable_equations;

  std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1>> solution;

  stokes() = default;
  stokes(femib::finite_element_space::finite_element_space<T, d, d> v,
         femib::finite_element_space::finite_element_space<T, d, 1> q,
         const femib::gauss::rule<T, d> &rule, T rho = 1.0, T mu = 1.0,
         std::function<femib::types::dvec<T, d>(femib::types::dvec<T, d>, T)>
             force = default_force<T, d>());

  std::vector<std::pair<types::dvec<T, d>, types::dvec<T, d>>>
  plot_velocity(T delta) {
    return V.plot(solution.back().topRows(V.size()), delta);
  }
  std::vector<std::pair<types::dvec<T, d>, types::dvec<T, 1>>>
  plot_pressure(T delta) {
    return Q.plot(solution.back().tail(Q.size()), delta);
  }
};

template <typename T, int d>
femib::util::sparse_solvable_equations<T>
augment_with_pressure_gauge(const stokes<T, d> &s) {

  femib::util::sparse_solvable_equations<T> base =
      femib::util::remove_edges<T>(s.AA, s.ff, s.bV, s.not_edges);

  int n = static_cast<int>(s.not_edges.size());
  Eigen::Matrix<T, Eigen::Dynamic, 1> constraint_row_reduced =
      Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(n);
  for (int k = 0; k < n; ++k) {
    int global_i = s.not_edges[k];
    if (global_i >= s.V.size()) {
      int pressure_i = global_i - s.V.size();
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

  // Tikhonov regularization
  T reg = T(1e-8) * s.A.diagonal().cwiseAbs().maxCoeff();
  if (reg == T(0)) {
    SPDLOG_WARN("femib::stokes_steady::augment_with_pressure_gauge: s.A's "
                "diagonal is all-zero: system may be left singular");
  }
  int nV_total = s.V.size();
  for (int k = 0; k < n; ++k) {
    bool is_velocity = s.not_edges[k] < nV_total;
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

template <typename T, int d> void rebuild_system(stokes<T, d> &s) {
  s.solvable_equations = augment_with_pressure_gauge<T, d>(s);
}

template <typename T, int d>
void rebuild_system(stokes<T, d> &s, const femib::gauss::rule<T, d> &rule) {
  std::function<T(femib::types::dvec<T, d>)> zero_boundary =
      [](const femib::types::dvec<T, d> &) { return T(0); };
  s.bV = femib::util::triplets2dense(
      femib::util::build_edges<T, d, d>(s.V, zero_boundary),
      s.V.size() + s.Q.size(), 1);
  s.pressure_constraint_vector =
      build_pressure_constraint_vector<T, d>(s.Q, rule);
  s.not_edges = femib::util::build_not_edges<T, d, d>(s.V);
  for (int i = 0; i < s.Q.size(); ++i) {
    s.not_edges.push_back(s.V.size() + i);
  }
  s.AA = assemble_saddle_point_matrix<T>(s.A, s.B, s.V.size(), s.Q.size());
  rebuild_system<T, d>(s);
}

template <typename T, int d, int e>
Eigen::Matrix<T, Eigen::Dynamic, 1> solve(const stokes<T, d> &s) {

  Eigen::BiCGSTAB<Eigen::SparseMatrix<T>, Eigen::IncompleteLUT<T>> solver;
  solver.setTolerance(1e-10);
  solver.setMaxIterations(500);
  solver.compute(s.solvable_equations.A);

  Eigen::Matrix<T, Eigen::Dynamic, 1> x = solver.solve(s.solvable_equations.b);
  if (solver.info() != Eigen::Success) {
    SPDLOG_ERROR("femib::stokes_steady::solve: BiCGSTAB+ILUT failed to "
                 "converge (Eigen::ComputationInfo = {}, iterations = {}, "
                 "estimated error = {})",
                 static_cast<int>(solver.info()), solver.iterations(),
                 solver.error());
  }

  // drop the last row, constraint on pressure
  Eigen::Matrix<T, Eigen::Dynamic, 1> x_no_lambda = x.topRows(x.rows() - 1);

  return femib::util::add_edges<T>(x_no_lambda, s.bV, s.V.size() + s.Q.size(),
                                   s.not_edges, s.V.nodes.E);
}

template <typename T, int d>
stokes<T, d>::stokes(
    femib::finite_element_space::finite_element_space<T, d, d> v,
    femib::finite_element_space::finite_element_space<T, d, 1> q,
    const femib::gauss::rule<T, d> &rule, T rho, T mu,
    std::function<femib::types::dvec<T, d>(femib::types::dvec<T, d>, T)> force)
    : V(std::move(v)), Q(std::move(q)), rho(rho), mu(mu),
      force(std::move(force)) {
  // the parameters rho, mu, force shadow the members (force is moved from)
  const T viscosity = this->mu;
  A = femib::util::triplets2sparse(
      femib::util::build_diagonal_matrix<T, d, d>(
          V, rule,
          [viscosity](femib::types::F<T, d, d> u, femib::types::F<T, d, d> w) {
            return stokes_a<T, d>(u, w, viscosity);
          }),
      V.size(), V.size());
  B = femib::util::triplets2sparse(
      femib::util::build_off_diagonal_matrix<T, d>(V, Q, rule, stokes_b<T, d>),
      V.size(), Q.size());

  std::function<femib::types::dvec<T, d>(femib::types::dvec<T, d>)> force_at_0 =
      [f = this->force](const femib::types::dvec<T, d> &x) {
        return f(x, T(0));
      };
  ff = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(V.size() + Q.size(), 1);
  ff.block(0, 0, V.size(), 1) = femib::util::build_load_vector<T, d, d>(
      V, rule, [force_at_0](femib::types::F<T, d, d> a) {
        return external_force<T, d>(a, force_at_0);
      });

  rebuild_system<T, d>(*this, rule);
}

} // namespace femib::stokes_steady
#endif
