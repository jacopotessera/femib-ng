#ifndef FEMIB_POISSON_HPP_INCLUDED_
#define FEMIB_POISSON_HPP_INCLUDED_

#include "../affine/affine.hpp"
#include "../finite_element_space/finite_element_space.hpp"
#include "../gauss/gauss.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <functional>
#include <iostream>
#include <vector>

#include "femib.hpp"

namespace femib::poisson {

// TODO still dense
template <typename T, int d, int e> struct poisson {

  femib::finite_element_space::finite_element_space<T, d, e> V;
  std::function<std::vector<Eigen::Triplet<T>>(
      femib::finite_element_space::finite_element_space<T, d, e> a,
      femib::finite_element_space::finite_element_space<T, d, e> b)>
      M;
  std::function<std::vector<Eigen::Triplet<T>>(
      femib::finite_element_space::finite_element_space<T, d, e> a)>
      f;

  Eigen::Matrix<T, Eigen::Dynamic, 1> dB;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> dM;
  Eigen::Matrix<T, Eigen::Dynamic, 1> dF;
};

// TODO add to constructor
template <typename T, int d, int e>
std::function<T(femib::types::dvec<T, d>)>
external_force(femib::types::F<T, d, e> a) {
  return [a](femib::types::dvec<T, d> x) {
    femib::types::dvec<T, e> value = a.x(x);
    T sum = 0;
    for (int k = 0; k < e; ++k) {
      sum += value(k);
    }
    return sum;
  };
}

template <typename T>
femib::util::solvable_equations<T>
remove_edges(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> dM,
             Eigen::Matrix<T, Eigen::Dynamic, 1> dF,
             Eigen::Matrix<T, Eigen::Dynamic, 1> dB, int rows,
             std::vector<int> not_edges) {

  Eigen::Matrix<T, Eigen::Dynamic, 1> ss =
      dM(Eigen::placeholders::all, Eigen::placeholders::all) *
      dB(Eigen::placeholders::all, Eigen::placeholders::all);

  Eigen::Matrix<T, Eigen::Dynamic, 1> bbb = (dF - ss)(not_edges, 0);
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> AAA =
      (dM)(not_edges, not_edges);

  return {AAA, bbb};
}

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> add_edges(

    Eigen::Matrix<T, Eigen::Dynamic, 1> xxx,
    Eigen::Matrix<T, Eigen::Dynamic, 1> dB, int rows,
    std::vector<int> not_edges, std::vector<int> nodesE) {

  Eigen::Matrix<T, Eigen::Dynamic, 1> xx;
  xx.resize(rows, 1);

  for (int i = 0; i < rows; i++) {
    xx(i, 0) = 0.0;
    auto k = std::find(not_edges.begin(), not_edges.end(), i);
    if (k != not_edges.end()) {
      xx(i, 0) = xxx(k - not_edges.begin(), 0);
    }
    auto kk = std::find(nodesE.begin(), nodesE.end(), i);
    if (kk != nodesE.end()) {
      xx(i, 0) = dB(i, 0);
    }
  }

  return xx;
}

// TODO constructor
template <typename T, int d, int e>
void init(poisson<T, d, e> &s, const femib::gauss::rule<T, d> &rule) {

  std::function<T(femib::types::dvec<T, d>)> b =
      [](const femib::types::dvec<T, d> &x) { return 0; };

  std::vector<Eigen::Triplet<T>> B = femib::util::build_edges<T, d, e>(s.V, b);

  std::vector<int> not_edges = femib::util::build_not_edges<T, d, e>(s.V);

  s.dB = femib::util::triplets2dense<T>(B, s.V.nodes.P.size(), 1);
  s.dM = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(
      util::build_stiffness_matrix(s.V, rule));
  s.dF = util::build_load_vector(s.V, rule, external_force<T, d, e>);
}

// TODO use BICSTAG
template <typename T, int d, int e>
Eigen::Matrix<T, Eigen::Dynamic, 1> solve(const poisson<T, d, e> &poisson) {

  std::vector<int> not_edges = femib::util::build_not_edges<T, d, e>(poisson.V);

  femib::util::solvable_equations<T> solvable_equations = remove_edges<T>(
      poisson.dM, poisson.dF, poisson.dB, poisson.V.nodes.P.size(), not_edges);

  Eigen::Matrix<T, Eigen::Dynamic, 1> x =
      solvable_equations.A.colPivHouseholderQr().solve(solvable_equations.b);

  return add_edges<T>(x, poisson.dB, poisson.V.nodes.P.size(), not_edges,
                      poisson.V.nodes.E);
}

} // namespace femib::poisson
#endif
