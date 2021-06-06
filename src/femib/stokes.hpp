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
stokes_a(const femib::types::F<T, d, d> &u, const femib::types::F<T, d, d> &v) {
  return [&u, &v](const femib::types::dvec<T, d> &x) { return dpi(u, v)(x); };
}

template <typename T, int d>
std::function<T(femib::types::dvec<T, d>)>
stokes_b(const femib::types::F<T, d, d> &u, const femib::types::F<T, d, 1> &q) {
  return [&u, &q](const femib::types::dvec<T, d> &x) {
    return div(u)(x) * q.x(x)(0);
  };
}

template <typename T, int d>
std::function<T(femib::types::dvec<T, d>)>
external_force(femib::types::F<T, d, d> a) {
  return
      [a](const femib::types::dvec<T, d> &x) { return a.x(x)[0] + a.x(x)[1]; };
}

template <typename T>
femib::util::solvable_equations<T>
remove_edges(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> dM,
             Eigen::Matrix<T, Eigen::Dynamic, 1> dF,
             Eigen::Matrix<T, Eigen::Dynamic, 1> bV,
             Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> bQ,

             std::vector<int> not_edges) {

  Eigen::Matrix<T, Eigen::Dynamic, 1> ss =
      dM(Eigen::all, Eigen::all) * bV(Eigen::all, Eigen::all);

  Eigen::Matrix<T, Eigen::Dynamic, 1> bbb = (dF - ss)(not_edges, 0);
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> AAA =
      (dM - bQ)(not_edges, not_edges);

  return {AAA, bbb};
}

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> add_edges(

    Eigen::Matrix<T, Eigen::Dynamic, 1> xxx,
    Eigen::Matrix<T, Eigen::Dynamic, 1> bV, int rowsV, int rowsQ,
    std::vector<int> not_edges, std::vector<int> nodesE,
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> bQ) {

  Eigen::Matrix<T, Eigen::Dynamic, 1> xx;
  xx.resize(rowsV + rowsQ, 1);

  for (int i = 0; i < rowsV; i++) {
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

  Eigen::Matrix<T, Eigen::Dynamic, 1> ppp = xxx.bottomRows(rowsQ - 1);

  Eigen::Matrix<T, 1, Eigen::Dynamic> PPP = bQ.block(0, 1, 1, rowsQ - 1);

  for (int i = rowsV; i < rowsV + rowsQ; i++) {
    xx(i, 0) = 0.0;
    if (i == rowsV) {
      xx(i, 0) = PPP * ppp;
    } else {
      auto k = std::find(not_edges.begin(), not_edges.end(), i);
      xx(i, 0) = xxx(k - not_edges.begin(), 0);
    }
  }

  return xx;
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

  std::function<T(femib::types::dvec<T, d>)> b =
      [](const femib::types::dvec<T, d> &x) { return 0.0; };

  s.bV =
      femib::util::triplets2dense(femib::util::build_edges<T, d, d>(s.V, b),
                                  s.V.nodes.P.size() + s.Q.nodes.P.size(), 1);

  s.bQ = femib::util::build_zero_mean_edges<T, d>(s.Q, rule);

  s.AA = Eigen::ArrayXXf::Zero(s.V.nodes.P.size() + s.Q.nodes.P.size(),
                               s.V.nodes.P.size() + s.Q.nodes.P.size());
  s.AA.block(0, 0, s.V.nodes.P.size(), s.V.nodes.P.size()) = s.A;

  s.AA.block(0, s.V.nodes.P.size(), s.V.nodes.P.size(), s.Q.nodes.P.size()) =
      s.B;

  s.AA.block(s.V.nodes.P.size(), 0, s.Q.nodes.P.size(), s.V.nodes.P.size()) =
      s.B.transpose();

  s.ff = Eigen::ArrayXXf::Zero(s.V.nodes.P.size() + s.Q.nodes.P.size(), 1);

  s.ff.block(0, 0, s.V.nodes.P.size(), 1) =
      femib::util::triplets2dense(result.F, s.V.nodes.P.size(), 1);

  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> bbQ =
      Eigen::ArrayXXf::Zero(s.V.nodes.P.size() + s.Q.nodes.P.size(),
                            s.V.nodes.P.size() + s.Q.nodes.P.size());
  bbQ.block(s.V.nodes.P.size(), s.V.nodes.P.size(), s.Q.nodes.P.size(),
            s.Q.nodes.P.size()) = s.bQ;

  std::vector<int> not_edges = femib::util::build_not_edges<T, d, d>(s.V);
  for (int i = 1; i < s.Q.nodes.P.size(); ++i) {
    not_edges.push_back(s.V.nodes.P.size() + i);
  }

  s.solvable_equations = remove_edges<T>(s.AA, s.ff, s.bV, bbQ, not_edges);
}

template <typename T, int d, int e>
Eigen::Matrix<T, Eigen::Dynamic, 1> solve(const stokes<T, d> &stokes) {

  Eigen::Matrix<T, Eigen::Dynamic, 1> x =
      stokes.solvable_equations.A.colPivHouseholderQr().solve(
          stokes.solvable_equations.b);

  std::vector<int> not_edges = femib::util::build_not_edges<T, d, d>(stokes.V);
  for (int i = 1; i < stokes.Q.nodes.P.size(); ++i) {
    not_edges.push_back(stokes.V.nodes.P.size() + i);
  }

  Eigen::Matrix<T, Eigen::Dynamic, 1> xx = add_edges<T>(
      x, stokes.bV, stokes.V.nodes.P.size(), stokes.Q.nodes.P.size(), not_edges,
      stokes.V.nodes.E, stokes.bQ);

  return xx;
}

} // namespace femib::stokes
#endif
