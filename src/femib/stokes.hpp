#ifndef FEMIB_STOKES_HPP_INCLUDED_
#define FEMIB_STOKES_HPP_INCLUDED_

#include "../femib/femib.hpp"
#include "../finite_element_space/finite_element_space.hpp"
#include "../gauss/gauss.hpp"
#include "../mesh/mesh.hpp"
#include "../types/differential_operation.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace femib::stokes {

template <typename T, int d> struct stokes {
  femib::finite_element_space::finite_element_space<T, d, d> V;
  femib::finite_element_space::finite_element_space<T, d, 1> Q;

  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> B;
  Eigen::Matrix<T, Eigen::Dynamic, 1> f;
  Eigen::Matrix<T, Eigen::Dynamic, 1> bV;

  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> bQ;

  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> AA;
  Eigen::Matrix<T, Eigen::Dynamic, 1> ff;

  femib::util::solvable_equations<T> solvable_equations;
};

template <typename T, int d>
std::function<T(femib::types::dvec<T, d>)>
stokes_a(femib::types::F<T, d, d> u, femib::types::F<T, d, d> v) {
  return [u, v](const femib::types::dvec<T, d> &x) {
    return T(2.0) * dpi(u, v)(x);
  }; // TODO  coefficient mu?
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
  return [a](const femib::types::dvec<T, d> &x) {
    return a.x(x)[0] + a.x(x)[1];
  }; // TODO same fix as poisson?
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
void rebuild_system(stokes<T, d> &s, const femib::gauss::rule<T, d> &rule) {

  std::function<T(femib::types::dvec<T, d>)> b =
      [](const femib::types::dvec<T, d> &x) { return 0.0; };

  s.bV =
      femib::util::triplets2dense(femib::util::build_edges<T, d, d>(s.V, b),
                                  s.V.nodes.P.size() + s.Q.nodes.P.size(), 1);

  Eigen::Matrix<T, 1, Eigen::Dynamic> domain_integral_row =
      femib::util::build_domain_integral_row<T, d>(s.Q, rule);

  s.AA = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(
      s.V.nodes.P.size() + s.Q.nodes.P.size(),
      s.V.nodes.P.size() + s.Q.nodes.P.size());
  s.AA.block(0, 0, s.V.nodes.P.size(), s.V.nodes.P.size()) = s.A;

  s.AA.block(0, s.V.nodes.P.size(), s.V.nodes.P.size(), s.Q.nodes.P.size()) =
      s.B;

  s.AA.block(s.V.nodes.P.size(), 0, s.Q.nodes.P.size(), s.V.nodes.P.size()) =
      s.B.transpose();

  std::vector<int> not_edges = femib::util::build_not_edges<T, d, d>(s.V);
  for (int i = 0; i < s.Q.nodes.P.size(); ++i) {
    not_edges.push_back(s.V.nodes.P.size() + i);
  }

  femib::util::solvable_equations<T> base =
      remove_edges<T>(s.AA, s.ff, s.bV, not_edges);

  // Augment with pressure gauge
  int n = not_edges.size();
  Eigen::Matrix<T, Eigen::Dynamic, 1> constraint_row_reduced =
      Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(n);
  for (int k = 0; k < n; ++k) {
    int global_i = not_edges[k];
    if (global_i >= s.V.nodes.P.size()) {
      int pressure_i = global_i - s.V.nodes.P.size();
      constraint_row_reduced(k) = domain_integral_row(pressure_i);
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

  s.solvable_equations = {AAA_aug, bbb_aug};
}

template <typename T, int d>
void init(stokes<T, d> &s, const femib::gauss::rule<T, d> &rule) {

  femib::util::build_diagonal_result<T> result =
      femib::util::build_diagonal<T, d, d>(s.V, rule, stokes_a<T, d>,
                                           external_force<T, d>);
  s.A = femib::util::triplets2dense(result.M, s.V.nodes.P.size(),
                                    s.V.nodes.P.size());
  s.B = femib::util::triplets2dense(
      femib::util::build_non_diagonal<T, d>(s.V, s.Q, rule, stokes_b<T, d>),
      s.V.nodes.P.size(), s.Q.nodes.P.size());

  s.ff = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(
      s.V.nodes.P.size() + s.Q.nodes.P.size(), 1);
  s.ff.block(0, 0, s.V.nodes.P.size(), 1) =
      femib::util::triplets2dense(result.F, s.V.nodes.P.size(), 1);

  rebuild_system<T, d>(s, rule);
}

template <typename T, int d, int e>
Eigen::Matrix<T, Eigen::Dynamic, 1> solve(const stokes<T, d> &stokes) {

  Eigen::Matrix<T, Eigen::Dynamic, 1> x =
      stokes.solvable_equations.A.colPivHouseholderQr().solve(
          stokes.solvable_equations.b); // TODO preconditioner?

  // drop the last row, constraint on pressure
  Eigen::Matrix<T, Eigen::Dynamic, 1> x_no_lambda = x.topRows(x.rows() - 1);

  std::vector<int> not_edges = femib::util::build_not_edges<T, d, d>(stokes.V);
  for (int i = 0; i < stokes.Q.nodes.P.size(); ++i) {
    not_edges.push_back(stokes.V.nodes.P.size() + i);
  }

  Eigen::Matrix<T, Eigen::Dynamic, 1> xx =
      add_edges<T>(x_no_lambda, stokes.bV, stokes.V.nodes.P.size(),
                   stokes.Q.nodes.P.size(), not_edges, stokes.V.nodes.E);

  return xx;
}

} // namespace femib::stokes
#endif
