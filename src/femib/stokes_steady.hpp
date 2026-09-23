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
#include <spdlog/spdlog.h>

namespace femib::stokes {

template <typename T, int d> struct stokes {
  femib::finite_element_space::finite_element_space<T, d, d> V;
  femib::finite_element_space::finite_element_space<T, d, 1> Q;
  T rho = 1.0;
  T mu = 1.0;

  Eigen::SparseMatrix<T> A;
  Eigen::SparseMatrix<T> B;
  Eigen::Matrix<T, Eigen::Dynamic, 1> f;
  Eigen::Matrix<T, Eigen::Dynamic, 1> bV;

  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> bQ;

  Eigen::SparseMatrix<T> AA;
  Eigen::Matrix<T, Eigen::Dynamic, 1> ff;

  femib::util::sparse_solvable_equations<T> solvable_equations;
};

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

template <typename T, int d>
std::function<T(femib::types::dvec<T, d>)>
external_force(femib::types::F<T, d, d> a) {
  return
      [a](const femib::types::dvec<T, d> &x) { return a.x(x)[0] + a.x(x)[1]; };
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

template <typename T, int d>
femib::util::sparse_solvable_equations<T> augment_with_pressure_gauge(
    const stokes<T, d> &s, const Eigen::SparseMatrix<T> &AA,
    const Eigen::Matrix<T, Eigen::Dynamic, 1> &ff,
    const Eigen::Matrix<T, Eigen::Dynamic, 1> &pressure_constraint_vector,
    const std::vector<int> &not_edges) {

  femib::util::sparse_solvable_equations<T> base =
      femib::util::remove_edges<T>(AA, ff, s.bV, not_edges);

  int n = static_cast<int>(not_edges.size());
  Eigen::Matrix<T, Eigen::Dynamic, 1> constraint_row_reduced =
      Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(n);
  for (int k = 0; k < n; ++k) {
    int global_i = not_edges[k];
    if (global_i >= s.V.nodes.P.size()) {
      int pressure_i = global_i - s.V.nodes.P.size();
      constraint_row_reduced(k) = pressure_constraint_vector(pressure_i);
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

  T reg = T(1e-8) * s.A.diagonal().cwiseAbs().maxCoeff();
  if (reg == T(0)) {
    SPDLOG_WARN("femib::stokes::augment_with_pressure_gauge: s.A's diagonal "
                "is all-zero: system may be left singular");
  }
  int nV_total = (int)s.V.nodes.P.size();
  for (int k = 0; k < n; ++k) {
    bool is_velocity = not_edges[k] < nV_total;
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

template <typename T, int d>
void rebuild_system(stokes<T, d> &s, const femib::gauss::rule<T, d> &rule) {

  std::function<T(femib::types::dvec<T, d>)> b =
      [](const femib::types::dvec<T, d> &x) { return 0.0; };

  s.bV =
      femib::util::triplets2dense(femib::util::build_edges<T, d, d>(s.V, b),
                                  s.V.nodes.P.size() + s.Q.nodes.P.size(), 1);

  Eigen::Matrix<T, Eigen::Dynamic, 1> pressure_constraint_vector =
      build_pressure_constraint_vector<T, d>(s.Q, rule);

  int nV = s.V.nodes.P.size();
  int nQ = s.Q.nodes.P.size();
  s.AA = assemble_saddle_point_matrix<T>(s.A, s.B, nV, nQ);

  std::vector<int> not_edges = femib::util::build_not_edges<T, d, d>(s.V);
  for (int i = 0; i < s.Q.nodes.P.size(); ++i) {
    not_edges.push_back(s.V.nodes.P.size() + i);
  }

  s.solvable_equations = augment_with_pressure_gauge<T, d>(
      s, s.AA, s.ff, pressure_constraint_vector, not_edges);
}

template <typename T, int d>
void init(stokes<T, d> &s, const femib::gauss::rule<T, d> &rule) {

  T mu = s.mu;
  s.A = femib::util::triplets2sparse(
      femib::util::build_diagonal_matrix<T, d, d>(
          s.V, rule,
          [mu](femib::types::F<T, d, d> u, femib::types::F<T, d, d> v) {
            return stokes_a<T, d>(u, v, mu);
          }),
      s.V.nodes.P.size(), s.V.nodes.P.size());
  s.B =
      femib::util::triplets2sparse(femib::util::build_off_diagonal_matrix<T, d>(
                                       s.V, s.Q, rule, stokes_b<T, d>),
                                   s.V.nodes.P.size(), s.Q.nodes.P.size());

  s.ff = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(
      s.V.nodes.P.size() + s.Q.nodes.P.size(), 1);
  s.ff.block(0, 0, s.V.nodes.P.size(), 1) =
      femib::util::build_load_vector<T, d, d>(s.V, rule, external_force<T, d>);

  rebuild_system<T, d>(s, rule);
}

template <typename T, int d, int e>
Eigen::Matrix<T, Eigen::Dynamic, 1> solve(const stokes<T, d> &stokes) {

  Eigen::BiCGSTAB<Eigen::SparseMatrix<T>, Eigen::IncompleteLUT<T>> solver;
  solver.setTolerance(1e-10);
  solver.setMaxIterations(500);
  solver.compute(stokes.solvable_equations.A);

  Eigen::Matrix<T, Eigen::Dynamic, 1> x =
      solver.solve(stokes.solvable_equations.b);
  if (solver.info() != Eigen::Success) {
    SPDLOG_ERROR("femib::stokes::solve: BiCGSTAB+ILUT failed to converge "
                 "(Eigen::ComputationInfo = {}, iterations = {}, "
                 "estimated error = {})",
                 static_cast<int>(solver.info()), solver.iterations(),
                 solver.error());
  }

  Eigen::Matrix<T, Eigen::Dynamic, 1> x_no_lambda = x.topRows(x.rows() - 1);

  std::vector<int> not_edges = femib::util::build_not_edges<T, d, d>(stokes.V);
  for (int i = 0; i < stokes.Q.nodes.P.size(); ++i) {
    not_edges.push_back(stokes.V.nodes.P.size() + i);
  }

  Eigen::Matrix<T, Eigen::Dynamic, 1> xx = femib::util::add_edges<T>(
      x_no_lambda, stokes.bV, stokes.V.nodes.P.size() + stokes.Q.nodes.P.size(),
      not_edges, stokes.V.nodes.E);

  return xx;
}

} // namespace femib::stokes
#endif
