#ifndef FEMIB_STOKES_HPP_INCLUDED_
#define FEMIB_STOKES_HPP_INCLUDED_

#include "../cuda/cuda.h"
#include "../femib/femib.hpp"
#include "../finite_element_space/finite_element_space.hpp"
#include "../gauss/gauss.hpp"
#include "../mesh/mesh.hpp"
#include "../types/differential_operation.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace femib::stokes_t {

template <typename T, int d> struct stokes {
  femib::finite_element_space::finite_element_space<T, d, d> V;
  femib::finite_element_space::finite_element_space<T, d, 1> Q;

  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> B;
  Eigen::Matrix<T, Eigen::Dynamic, 1> f;
  Eigen::Matrix<T, Eigen::Dynamic, 1> bV;
  femib::util::build_diagonal_result<T> result;

  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> bQ;

  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> AA;
  Eigen::Matrix<T, Eigen::Dynamic, 1> ff;

  femib::util::solvable_equations<T> solvable_equations;

  T deltat = 0.1;
  std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1>> solution;
  std::vector<std::vector<std::vector<std::vector<float>>>> plot;
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
      xx(i, 0) = -PPP * ppp;
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

  s.result = result;
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
      femib::util::triplets2dense(s.result.F, s.V.nodes.P.size(), 1);

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

template <typename T, int d> void advance(stokes<T, d> &s) {
  Eigen::Matrix<T, Eigen::Dynamic, 1> u_1;
  if (s.solution.size() == 0)
    u_1 = Eigen::ArrayXf::Zero(s.V.nodes.P.size(), 1);
  else
    u_1 = s.solution[s.solution.size() - 1];
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> DD =
      (1 / s.deltat) *
      Eigen::MatrixXf::Identity(s.V.nodes.P.size(), s.V.nodes.P.size());
  Eigen::Matrix<T, Eigen::Dynamic, 1> dd = (1 / s.deltat) * u_1;

  s.AA.block(0, 0, s.V.nodes.P.size(), s.V.nodes.P.size()) = s.A + DD;

  s.ff.block(0, 0, s.V.nodes.P.size(), 1) =
      femib::util::triplets2dense(s.result.F, s.V.nodes.P.size(), 1) + dd;

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

  // TODO
  Eigen::Matrix<T, Eigen::Dynamic, 1> xx = solve<T, d, 1>(s);

  std::vector<std::vector<std::vector<float>>> uuu;

  femib::types::mesh mesh = s.V.mesh;
  femib::types::box<float, 2> box = femib::mesh::find_box<float, 2>(mesh);

  femib::types::box<float, 2> boxx =
      femib::mesh::lin_spaced<float, 2>(box, 0.027);
  bool N[boxx.size() * mesh.N.size()];

  femib::cuda::serial_accurate<float, 2>(boxx.data(), boxx.size(),
                                         mesh.N.data(), mesh.N.size(), N);

  std::vector<int> NNN;

  for (int i = 0; i < boxx.size(); ++i) {
    for (int n = 0; n < mesh.N.size(); ++n) {
      if (N[i * mesh.N.size() + n]) {
        NNN.push_back(n);
        break;
      }
    }
  }

  for (int i = 0; i < boxx.size(); ++i) {
    femib::types::dvec<T, d> point = boxx[i];
    femib::types::dvec<T, d> res = {0, 0};
    // if (NNN[i] < s.V.mesh.T.size()) {
    femib::types::dtrian<T, d> t = mesh.N[NNN[i]];

    for (int j = 0; j < s.V.finite_element.base_functions.size(); ++j) {

      femib::types::F<T, d, d> f = s.V.finite_element.base_functions[j];

      std::function<femib::types::dvec<T, d>(femib::types::dvec<T, d>)> g =
          [&](const femib::types::dvec<T, d> &x) {
            return xx(s.V.nodes.get_index(j, NNN[i])) *
                   f.x(femib::affine::affine_inv(t, x));
          };
      res += g(point);
    }
    //}

    std::vector<std::vector<float>> uuuu = {{point(0), point(1)},
                                            {res(0), res(1)}};
    uuu.emplace_back(uuuu);
  }

  s.plot.emplace_back(uuu);
  s.solution.emplace_back(xx);
}

template <typename T, int d, int e>
Eigen::Matrix<T, Eigen::Dynamic, 1> solve(const stokes<T, d> &s) {

  Eigen::Matrix<T, Eigen::Dynamic, 1> x =
      s.solvable_equations.A.colPivHouseholderQr().solve(
          s.solvable_equations.b);

  std::vector<int> not_edges = femib::util::build_not_edges<T, d, d>(s.V);
  for (int i = 1; i < s.Q.nodes.P.size(); ++i) {
    not_edges.push_back(s.V.nodes.P.size() + i);
  }

  Eigen::Matrix<T, Eigen::Dynamic, 1> xx =
      add_edges<T>(x, s.bV, s.V.nodes.P.size(), s.Q.nodes.P.size(), not_edges,
                   s.V.nodes.E, s.bQ);

  return xx;
}

} // namespace femib::stokes_t
#endif
