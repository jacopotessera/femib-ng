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

  Eigen::Matrix<T, 1, Eigen::Dynamic> bQ; // = 0

  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> AA;
  Eigen::Matrix<T, Eigen::Dynamic, 1> ff;
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
  return [a](femib::types::dvec<T, d> x) { return -1; };
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

  s.bV = femib::util::triplets2dense(femib::util::build_edges<T, d, d>(s.V, b),
                                     s.V.nodes.P.size(), 1);

  s.bQ = femib::util::triplets2dense(
      femib::util::build_zero_mean_edges<T, d>(s.Q, rule), 1,
      s.Q.nodes.P.size());

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
}

} // namespace femib::stokes
#endif
