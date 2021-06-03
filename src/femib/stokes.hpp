#ifndef FEMIB_STOKES_HPP_INCLUDED_
#define FEMIB_STOKES_HPP_INCLUDED_

#include "../finite_element_space/finite_element_space.hpp"
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
}

// template <typename T, int d>
// std::function<T(femib::types::dvec<T, d>)> a(femib::types::F<T, d, d> u,
//                                               femib::types::F<T, d, d> v) {
//  return [u, v](femib::types::dvec<T, d> x) {
//    return 0.5*( u.dx(x)+u.dx(x).transpose())* (v.dx(x)+v.dx(x).transpose());
//  };
//}

template <typename T, int d>
std::function<T(femib::types::dvec<T, d>)>
stokes_b(femib::types::F<T, d, d> u, femib::types::F<T, d, 1> q) {
  return [u, q](femib::types::dvec<T, d> x) { return div(u) * q; };
}

template <typename T, int d>
std::vector<Eigen::Triplet<T>>
build_diagonal(femib::finite_element_space::finite_element_space<T, d, d> v,
               femib::finite_element_space::finite_element_space<T, d, 1> q,
               femib::gauss::rule<T, d> rule) {
  std::vector<Eigen::Triplet<T>> BB;
  for (int n = 0; n < v.mesh.T.size(); ++n) {
    femib::types::dtrian<T, d> t = s.mesh[n];
    for (int i = 0; i < v.finite_element.base_functions.size(); ++i) {
      femib::types::F<T, d, d> a;
      a.x =
          [&](const femib::types::dvec<T, d> &x) {
            return (v.finite_element.base_functions[i].x(
                femib::affine::affineBinv(t) *
                (x - femib::affine::affineb(t))));
          },
      a.dx = [&](const femib::types::dvec<T, d> &x) {
        return (femib::affine::affineBinv(t) *
                v.finite_element.base_functions[i].dx(
                    femib::affine::affineBinv(t) *
                    (x - femib::affine::affineb(t))));
      };
      for (int j = 0; j < q.finite_element.base_functions.size(); ++j) {
        femib::types::F<T, d, 1> b;
        b.x = [&](const femib::types::dvec<T, d> &x) {
          return (q.finite_element.base_functions[i].x(
              femib::affine::affineBinv(t) * (x - femib::affine::affineb(t))));
        };
        T m = femib::mesh::integrate<T, d>(rule, stokes_b(a, b), t);
        BB.push_back(
            Eigen::Triplet<T>(femib::types::get_index<T, d>(s.nodes, i, n),
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
void init(stokes<T, d, e> &s, const femib::gauss::rule<T, d> &rule) {}
        B = femib::util::triplets2dense<T>(build_diagonal(s.V, s.Q, rule);
        } // namespace femib::stokes
#endif
