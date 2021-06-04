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

// template <typename T, int d>
// std::function<T(femib::types::dvec<T, d>)> a(femib::types::F<T, d, d> u,
//                                               femib::types::F<T, d, d> v) {
//  return [u, v](femib::types::dvec<T, d> x) {
//    return 0.5*( u.dx(x)+u.dx(x).transpose())* (v.dx(x)+v.dx(x).transpose());
//  };
//}

template <typename T, int d>
std::function<T(femib::types::dvec<T, d>)>
stokes_b(const femib::types::F<T, d, d> &u, const femib::types::F<T, d, 1> &q) {
  // return project<T, d>(q.x, 0) * div(u);
  return [&u, &q](const femib::types::dvec<T, d> &x) {
    return div(u)(x) * q.x(x)(0);
  };
}

template <typename T, int d>
std::vector<Eigen::Triplet<T>> build_diagonal(
    const femib::finite_element_space::finite_element_space<T, d, d> &v,
    const femib::finite_element_space::finite_element_space<T, d, 1> &q,
    const femib::gauss::rule<T, d> &rule) {
  std::vector<Eigen::Triplet<T>> BB;
  for (int n = 0; n < v.mesh.T.size(); ++n) {
    femib::types::dtrian<T, d> t = v.mesh[n];
    for (int i = 0; i < v.finite_element.base_functions.size(); ++i) {
      femib::types::F<T, d, d> a =
          femib::util::base_function2real_function<T, d, d>(v, n, i);
      for (int j = 0; j < q.finite_element.base_functions.size(); ++j) {
        femib::types::F<T, d, 1> b;
        b.x = [&](const femib::types::dvec<T, d> &x) {
          return (q.finite_element.base_functions[j].x(
              femib::affine::affineBinv(t) * (x - femib::affine::affineb(t))));
        };
        T m = femib::mesh::integrate<T, d>(rule, stokes_b(a, b), t);
        BB.push_back(
            Eigen::Triplet<T>(femib::types::get_index<T, d>(v.nodes, i, n),
                              femib::types::get_index<T, d>(q.nodes, j, n), m));
      }
      // T f_ = femib::mesh::integrate<T, d>(rule, ggg(a), t);
      // FF.push_back(Eigen::Triplet<T>(
      //    femib::types::get_index<T, d>(s.nodes, i, n), 0, f_));
    }
  }
  return BB;
}

template <typename T, int d>
void init(stokes<T, d> &s, const femib::gauss::rule<T, d> &rule) {
  s.B = femib::util::triplets2dense(build_diagonal(s.V, s.Q, rule),
                                    s.V.nodes.P.size(), s.Q.nodes.P.size());
}

} // namespace femib::stokes
#endif
