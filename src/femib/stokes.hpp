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
  Eigen::Matrix<T, Eigen::Dynamic, 1> b;
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
  return [a](femib::types::dvec<T, d> x) {
    return -400 * x(0) * x(1) * a.x(x)[0] + 10 * a.x(x)[1];
  };
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
}

} // namespace femib::stokes
#endif
