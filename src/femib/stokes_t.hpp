#ifndef FEMIB_STOKES_T_HPP_INCLUDED_
#define FEMIB_STOKES_T_HPP_INCLUDED_

#include "../femib/femib.hpp"
#include "../finite_element_space/finite_element_space.hpp"
#include "../gauss/gauss.hpp"
#include "../mesh/mesh.hpp"
#include "../types/differential_operation.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <spdlog/spdlog.h>

namespace femib::stokes_t {

// TODO constructor, V, Q, deltat, gauss-rule, force, etc
template <typename T, int d> struct stokes {
  femib::finite_element_space::finite_element_space<T, d, d> V;
  femib::finite_element_space::finite_element_space<T, d, 1> Q;
  femib::gauss::rule<T, d> rule;
  std::function<femib::types::dvec<T, d>(femib::types::dvec<T, d>, T)> force =
      [](femib::types::dvec<T, d>, T) -> Eigen::Matrix<float, 2, 1> {
    return Eigen::Matrix<T, d, 1>::Zero();
  }; // external force f(x,t)
  T deltat = 0.1; // TODO eh

  Eigen::SparseMatrix<T> A;
  Eigen::SparseMatrix<T> M;
  Eigen::SparseMatrix<T> B;
  Eigen::Matrix<T, Eigen::Dynamic, 1> f;
  Eigen::Matrix<T, Eigen::Dynamic, 1> bV;
  femib::util::build_diagonal_result<T> result;

  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> bQ;
  Eigen::Matrix<T, 1, Eigen::Dynamic> domain_integral_row;

  Eigen::SparseMatrix<T> AA;
  Eigen::Matrix<T, Eigen::Dynamic, 1> ff;

  util::sparse_solvable_equations<T> solvable_equations;

  std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1>> solution;
  std::vector<std::vector<std::pair<types::dvec<T, d>, types::dvec<T, d>>>>
      plotV; // TODO we need this?? cant we calculate from solution if needed?
  std::vector<std::vector<std::pair<types::dvec<T, d>, types::dvec<T, 1>>>>
      plotQ; // TODO we need this?? cant we calculate from solution if needed?

  T time = 0;
};

// TODO unify with stokes
template <typename T, int d>
std::function<T(femib::types::dvec<T, d>)>
stokes_a(femib::types::F<T, d, d> u, femib::types::F<T, d, d> v) {
  return [u, v](const femib::types::dvec<T, d> &x) {
    return T(2.0) * dpi(u, v)(x);
  }; // TODO coefficient mu?
}

template <typename T, int d>
std::function<T(femib::types::dvec<T, d>)>
stokes_b(femib::types::F<T, d, d> u, femib::types::F<T, d, 1> q) {
  return [u, q](const femib::types::dvec<T, d> &x) {
    return T(-1.0) * div(u)(x) * q.x(x)(0);
  };
}

template <typename T, int d>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
mass_matrix(const femib::finite_element_space::finite_element_space<T, d, d> &V,
            const femib::gauss::rule<T, d> &rule) {
  femib::util::build_diagonal_result<T> result =
      femib::util::build_diagonal<T, d, d>(V, rule, mass<T, d>, zero<T, d>);
  return femib::util::triplets2dense(result.M, V.nodes.P.size(),
                                     V.nodes.P.size());
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

// select not edges cols/rows
template <typename T>
Eigen::SparseMatrix<T> selection_matrix(const std::vector<int> &rows, int n) {
  Eigen::SparseMatrix<T> P(rows.size(), n);
  std::vector<Eigen::Triplet<T>> triplets;
  triplets.reserve(rows.size());
  for (size_t i = 0; i < rows.size(); ++i) {
    triplets.push_back(Eigen::Triplet<T>((int)i, rows[i], T(1)));
  }
  P.setFromTriplets(triplets.begin(), triplets.end());
  return P;
}

template <typename T>
util::sparse_solvable_equations<T>
remove_edges(const Eigen::SparseMatrix<T> &dM,
             const Eigen::Matrix<T, Eigen::Dynamic, 1> &dF,
             const Eigen::Matrix<T, Eigen::Dynamic, 1> &bV,
             const std::vector<int> &not_edges) {

  Eigen::Matrix<T, Eigen::Dynamic, 1> ss = dM * bV;

  int n = static_cast<int>(not_edges.size());
  Eigen::Matrix<T, Eigen::Dynamic, 1> bbb(n);
  for (int i = 0; i < n; ++i) {
    bbb(i) = dF(not_edges[i]) - ss(not_edges[i]);
  }

  Eigen::SparseMatrix<T> P =
      selection_matrix<T>(not_edges, static_cast<int>(dM.rows()));
  Eigen::SparseMatrix<T> AAA = P * dM * P.transpose();

  return {AAA, bbb};
}

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> add_edges(

    Eigen::Matrix<T, Eigen::Dynamic, 1> xxx,
    Eigen::Matrix<T, Eigen::Dynamic, 1> bV, int rowsV, int rowsQ,
    std::vector<int> not_edges, std::vector<int> nodesE) {

  Eigen::Matrix<T, Eigen::Dynamic, 1> xx;
  xx.resize(rowsV + rowsQ, 1);

  for (int i = 0; i < rowsV + rowsQ; i++) {
    xx(i, 0) = 0.0;
    auto k = std::find(not_edges.begin(), not_edges.end(), i);
    if (k != not_edges.end()) {
      xx(i, 0) = xxx(k - not_edges.begin(), 0);
    }
    auto kk = std::find(nodesE.begin(), nodesE.end(), i);
    if (kk != nodesE.end()) {
      xx(i, 0) = bV(i, 0);
    }
  }

  return xx;
}

template <typename T, int d>
std::vector<int> build_stokes_t_not_edges(const stokes<T, d> &s) {
  std::vector<int> not_edges = femib::util::build_not_edges<T, d, d>(s.V);
  for (int i = 0; i < s.Q.nodes.P.size(); ++i) {
    not_edges.push_back(s.V.nodes.P.size() + i);
  }
  return not_edges;
}

template <typename T, int d>
util::sparse_solvable_equations<T>
augment_with_pressure_gauge(const stokes<T, d> &s,
                            const Eigen::SparseMatrix<T> &AA,
                            const Eigen::Matrix<T, Eigen::Dynamic, 1> &ff,
                            const std::vector<int> &not_edges) {

  util::sparse_solvable_equations<T> base =
      remove_edges<T>(AA, ff, s.bV, not_edges);

  int n = static_cast<int>(not_edges.size());
  Eigen::Matrix<T, Eigen::Dynamic, 1> constraint_row_reduced =
      Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(n);
  for (int k = 0; k < n; ++k) {
    int global_i = not_edges[k];
    if (global_i >= s.V.nodes.P.size()) {
      int pressure_i = global_i - s.V.nodes.P.size();
      constraint_row_reduced(k) = s.domain_integral_row(pressure_i);
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
    SPDLOG_WARN(
        "femib::stokes_t::augment_with_pressure_gauge: s.A's diagonal is "
        "all-zero: system may be left singular");
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

// assembles the saddle-point matrix [[top_left, B], [B^T, 0]]
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

template <typename T, int d>
void init(stokes<T, d> &s, const femib::gauss::rule<T, d> &rule) {
  s.rule = rule;

  std::function<femib::types::dvec<T, d>(femib::types::dvec<T, d>)> force0 =
      [force = s.force, time0 = s.time](const femib::types::dvec<T, d> &x) {
        return force(x, time0);
      };
  auto ggg = [force0](femib::types::F<T, d, d> a) {
    return external_force<T, d>(a, force0);
  };

  femib::util::build_diagonal_result<T> result =
      femib::util::build_diagonal<T, d, d>(s.V, rule, stokes_a<T, d>, ggg);

  s.result = result;
  s.A = femib::util::triplets2sparse(result.M, s.V.nodes.P.size(),
                                     s.V.nodes.P.size());
  // TODO build sparse mass_matrix
  s.M = mass_matrix<T, d>(s.V, rule).sparseView();
  s.B = femib::util::triplets2sparse(
      femib::util::build_non_diagonal<T, d>(s.V, s.Q, rule, stokes_b<T, d>),
      s.V.nodes.P.size(), s.Q.nodes.P.size());

  std::function<T(femib::types::dvec<T, d>)> b =
      [](const femib::types::dvec<T, d> &x) { return 0.0; };

  s.bV =
      femib::util::triplets2dense(femib::util::build_edges<T, d, d>(s.V, b),
                                  s.V.nodes.P.size() + s.Q.nodes.P.size(), 1);

  s.domain_integral_row =
      femib::util::build_domain_integral_row<T, d>(s.Q, rule);

  int nV = s.V.nodes.P.size();
  int nQ = s.Q.nodes.P.size();
  s.AA = assemble_saddle_point_matrix<T>(s.A, s.B, nV, nQ);

  s.ff = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(nV + nQ, 1);

  s.ff.block(0, 0, nV, 1) = femib::util::triplets2dense(s.result.F, nV, 1);

  std::vector<int> not_edges = build_stokes_t_not_edges<T, d>(s);

  s.solvable_equations =
      augment_with_pressure_gauge<T, d>(s, s.AA, s.ff, not_edges);
}

// TODO we need to give better names to stuff...
// Rebuilds:
// * the time-dependent load vector (s.force evaluated at s.time,
// which the caller must have already advanced to t^{n+1});
// * dd (the backward-Euler mass term's known part);
// * extra_velocity_rhs.
// Then
// * re-derives the solvable system from the current s.AA/s.ff,
// extracted so the Picard loop for time-dependent Navier-Stokes
// can call it too.
template <typename T, int d>
void rebuild_system(
    stokes<T, d> &s, const Eigen::Matrix<T, Eigen::Dynamic, 1> &dd,
    std::optional<Eigen::Matrix<T, Eigen::Dynamic, 1>> extra_velocity_rhs) {

  std::function<femib::types::dvec<T, d>(femib::types::dvec<T, d>)> force_n1 =
      [force = s.force, time_n1 = s.time](const femib::types::dvec<T, d> &x) {
        return force(x, time_n1);
      };
  auto ggg = [force_n1](femib::types::F<T, d, d> a) {
    return external_force<T, d>(a, force_n1);
  };
  std::vector<Eigen::Triplet<T>> F_triplets =
      femib::util::build_vector<T, d, d>(s.V, s.rule, ggg);

  Eigen::Matrix<T, Eigen::Dynamic, 1> velocity_rhs =
      femib::util::triplets2dense(F_triplets, s.V.nodes.P.size(), 1) + dd;

  if (extra_velocity_rhs.has_value()) {
    velocity_rhs += extra_velocity_rhs.value();
  }

  s.ff.block(0, 0, s.V.nodes.P.size(), 1) = velocity_rhs;

  // The velocity-velocity block of AA/ff just changed (the backward-Euler
  // mass term, and possibly a Picard-updated convection term), so the
  // solvable system must be rebuilt.
  // s.domain_integral_row itself does not change, so it is reused.
  std::vector<int> not_edges = build_stokes_t_not_edges<T, d>(s);

  s.solvable_equations =
      augment_with_pressure_gauge<T, d>(s, s.AA, s.ff, not_edges);
}

template <typename T, int d, int e>
Eigen::Matrix<T, Eigen::Dynamic, 1> solve(const stokes<T, d> &s) {

  Eigen::SparseLU<Eigen::SparseMatrix<T>> solver;
  solver.compute(s.solvable_equations.A);
  if (solver.info() != Eigen::Success) {
    SPDLOG_ERROR("femib::stokes_t::solve: SparseLU factorization failed "
                 "(Eigen::ComputationInfo = {})",
                 static_cast<int>(solver.info()));
  }

  Eigen::Matrix<T, Eigen::Dynamic, 1> x = solver.solve(s.solvable_equations.b);

  // drop the last row, constraint on pressure
  Eigen::Matrix<T, Eigen::Dynamic, 1> x_no_lambda = x.topRows(x.rows() - 1);

  std::vector<int> not_edges = build_stokes_t_not_edges<T, d>(s);

  Eigen::Matrix<T, Eigen::Dynamic, 1> xx =
      add_edges<T>(x_no_lambda, s.bV, s.V.nodes.P.size(), s.Q.nodes.P.size(),
                   not_edges, s.V.nodes.E);

  return xx;
}

template <typename T, int d>
void advance(stokes<T, d> &s, std::optional<Eigen::Matrix<T, Eigen::Dynamic, 1>>
                                  extra_velocity_rhs = std::nullopt) {
  s.time += s.deltat;

  Eigen::Matrix<T, Eigen::Dynamic, 1> u_1;
  if (s.solution.size() == 0)
    u_1 = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(s.V.nodes.P.size(), 1);
  else
    // last timestep velocity
    u_1 = s.solution[s.solution.size() - 1].topRows(s.V.nodes.P.size());

  Eigen::SparseMatrix<T> DD = (1 / s.deltat) * s.M;
  Eigen::Matrix<T, Eigen::Dynamic, 1> dd = (1 / s.deltat) * (s.M * u_1);

  int nV = s.V.nodes.P.size();
  int nQ = s.Q.nodes.P.size();
  Eigen::SparseMatrix<T> top_left = s.A + DD;
  s.AA = assemble_saddle_point_matrix<T>(top_left, s.B, nV, nQ);

  rebuild_system<T, d>(s, dd, extra_velocity_rhs);

  Eigen::Matrix<T, Eigen::Dynamic, 1> xx = femib::stokes_t::solve<T, d, 1>(s);

  s.plotV.emplace_back(s.V.plot(xx.head(s.V.nodes.P.size()), 0.01)); // TODO
  s.plotQ.emplace_back(s.Q.plot(xx.tail(s.Q.nodes.P.size()), 0.01)); // TODO
  s.solution.emplace_back(xx);
}

} // namespace femib::stokes_t
#endif
