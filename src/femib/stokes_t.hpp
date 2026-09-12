#ifndef FEMIB_STOKES_T_HPP_INCLUDED_
#define FEMIB_STOKES_T_HPP_INCLUDED_

#include "../femib/femib.hpp"
#include "../finite_element_space/finite_element_space.hpp"
#include "../gauss/gauss.hpp"
#include "../mesh/mesh.hpp"
#include "../types/differential_operation.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>

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

  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> M;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> B;
  Eigen::Matrix<T, Eigen::Dynamic, 1> f;
  Eigen::Matrix<T, Eigen::Dynamic, 1> bV;
  femib::util::build_diagonal_result<T> result;

  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> bQ;
  Eigen::Matrix<T, 1, Eigen::Dynamic> domain_integral_row;

  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> AA;
  Eigen::Matrix<T, Eigen::Dynamic, 1> ff;

  femib::util::solvable_equations<T> solvable_equations;

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
mass_matrix(femib::finite_element_space::finite_element_space<T, d, d> &V,
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

template <typename T>
femib::util::solvable_equations<T>
remove_edges(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> dM,
             Eigen::Matrix<T, Eigen::Dynamic, 1> dF,
             Eigen::Matrix<T, Eigen::Dynamic, 1> bV,

             std::vector<int> not_edges) {

  Eigen::Matrix<T, Eigen::Dynamic, 1> ss =
      dM(Eigen::placeholders::all, Eigen::placeholders::all) *
      bV(Eigen::placeholders::all, Eigen::placeholders::all);

  Eigen::Matrix<T, Eigen::Dynamic, 1> bbb = (dF - ss)(not_edges, 0);
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> AAA =
      dM(not_edges, not_edges);

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
femib::util::solvable_equations<T> augment_with_pressure_gauge(
    const stokes<T, d> &s,
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &AA,
    const Eigen::Matrix<T, Eigen::Dynamic, 1> &ff,
    const std::vector<int> &not_edges) {

  femib::util::solvable_equations<T> base =
      remove_edges<T>(AA, ff, s.bV, not_edges);

  int n = not_edges.size();
  Eigen::Matrix<T, Eigen::Dynamic, 1> constraint_row_reduced =
      Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(n);
  for (int k = 0; k < n; ++k) {
    int global_i = not_edges[k];
    if (global_i >= s.V.nodes.P.size()) {
      int pressure_i = global_i - s.V.nodes.P.size();
      constraint_row_reduced(k) = s.domain_integral_row(pressure_i);
    }
  }

  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> AAA_aug =
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(n + 1, n + 1);
  AAA_aug.block(0, 0, n, n) = base.A;
  AAA_aug.block(0, n, n, 1) = constraint_row_reduced;
  AAA_aug.block(n, 0, 1, n) = constraint_row_reduced.transpose();

  Eigen::Matrix<T, Eigen::Dynamic, 1> bbb_aug =
      Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(n + 1);
  bbb_aug.topRows(n) = base.b;

  return {AAA_aug, bbb_aug};
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
  s.A = femib::util::triplets2dense(result.M, s.V.nodes.P.size(),
                                    s.V.nodes.P.size());
  s.M = mass_matrix<T, d>(s.V, rule);
  s.B = femib::util::triplets2dense(
      femib::util::build_non_diagonal<T, d>(s.V, s.Q, rule, stokes_b<T, d>),
      s.V.nodes.P.size(), s.Q.nodes.P.size());

  std::function<T(femib::types::dvec<T, d>)> b =
      [](const femib::types::dvec<T, d> &x) { return 0.0; };

  s.bV =
      femib::util::triplets2dense(femib::util::build_edges<T, d, d>(s.V, b),
                                  s.V.nodes.P.size() + s.Q.nodes.P.size(), 1);

  s.domain_integral_row =
      femib::util::build_domain_integral_row<T, d>(s.Q, rule);

  s.AA = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(
      s.V.nodes.P.size() + s.Q.nodes.P.size(),
      s.V.nodes.P.size() + s.Q.nodes.P.size());

  s.AA.block(0, 0, s.V.nodes.P.size(), s.V.nodes.P.size()) = s.A;

  s.AA.block(0, s.V.nodes.P.size(), s.V.nodes.P.size(), s.Q.nodes.P.size()) =
      s.B;

  s.AA.block(s.V.nodes.P.size(), 0, s.Q.nodes.P.size(), s.V.nodes.P.size()) =
      s.B.transpose();

  s.ff = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(
      s.V.nodes.P.size() + s.Q.nodes.P.size(), 1);

  s.ff.block(0, 0, s.V.nodes.P.size(), 1) =
      femib::util::triplets2dense(s.result.F, s.V.nodes.P.size(), 1);

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

  Eigen::Matrix<T, Eigen::Dynamic, 1> x =
      s.solvable_equations.A.colPivHouseholderQr().solve(
          s.solvable_equations.b);

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

  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> DD = (1 / s.deltat) * s.M;
  Eigen::Matrix<T, Eigen::Dynamic, 1> dd = (1 / s.deltat) * (s.M * u_1);

  s.AA.block(0, 0, s.V.nodes.P.size(), s.V.nodes.P.size()) = s.A + DD;

  rebuild_system<T, d>(s, dd, extra_velocity_rhs);

  Eigen::Matrix<T, Eigen::Dynamic, 1> xx = femib::stokes_t::solve<T, d, 1>(s);

  s.plotV.emplace_back(s.V.plot(xx.head(s.V.nodes.P.size()), 0.01)); // TODO
  s.plotQ.emplace_back(s.Q.plot(xx.tail(s.Q.nodes.P.size()), 0.01)); // TODO
  s.solution.emplace_back(xx);
}

} // namespace femib::stokes_t
#endif
