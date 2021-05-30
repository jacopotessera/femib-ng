#ifndef FEMIB_HPP_INCLUDED_
#define FEMIB_HPP_INCLUDED_

#include "../affine/affine.hpp"
#include "../finite_element_space/finite_element_space.hpp"
#include "../gauss/gauss.hpp"
#include "../mesh/mesh.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <functional>
#include <vector>

int get_index(const femib::types::nodes<float, 2> &nodes, int i, int n) {
  return nodes.T[n][i];
}

namespace femib::poisson {

template <typename T, int d, int e> struct poisson {

  femib::finite_element_space::finite_element_space<T, d, e> V;
  std::function<std::vector<Eigen::Triplet<T>>(
      femib::finite_element_space::finite_element_space<T, d, e> a,
      femib::finite_element_space::finite_element_space<T, d, e> b)>
      M;
  std::function<std::vector<Eigen::Triplet<T>>(
      femib::finite_element_space::finite_element_space<T, d, e> a)>
      f;
};

template <typename T, int d, int e>
std::function<T(femib::types::dvec<T, d>)>
ddot(femib::types::F<float, 2, 1> a, femib::types::F<float, 2, 1> b) {
  return [a, b](femib::types::dvec<float, 2> x) {
    return a.dx(x)[0] * b.dx(x)[0] + a.dx(x)[1] * b.dx(x)[1];
  };
}

template <typename T, int d, int e>
std::vector<Eigen::Triplet<T>>
build_diagonal(femib::finite_element_space::finite_element_space<T, d, e> s,
               femib::gauss::rule<T, d> rule,

               std::function<

                   std::function<T(femib::types::dvec<T, d>)>

                   (femib::types::F<float, 2, 1>, femib::types::F<float, 2, 1>)

                   >
                   fff) {

  std::vector<Eigen::Triplet<float>> MM;
  for (int n = 0; n < s.mesh.T.size(); ++n) {
    for (int i = 0; i < s.finite_element.base_functions.size(); ++i) {
      for (int j = 0; j < s.finite_element.base_functions.size(); ++j) {
        femib::types::F<float, 2, 1> a;
        femib::types::F<float, 2, 1> b;
        femib::types::dtrian<float, 2> t = s.mesh[n];
        a.dx = [&](const femib::types::dvec<float, 2> &x) {
          return (femib::affine::affineBinv(t) *
                  s.finite_element.base_functions[i].dx(
                      femib::affine::affineBinv(t) *
                      (x - femib::affine::affineb(t))));
        };
        b.dx = [&](const femib::types::dvec<float, 2> &x) {
          return (femib::affine::affineBinv(t) *
                  s.finite_element.base_functions[j].dx(
                      femib::affine::affineBinv(t) *
                      (x - femib::affine::affineb(t))));
        };
        float m = femib::mesh::integrate<float, 2>(
            rule,
            //[&t, &f, i, j, &a, &b](femib::types::dvec<float, 2> x) {
            //  return a.dx(x)[0] * b.dx(x)[0] + a.dx(x)[1] * b.dx(x)[1];
            //},
            fff(a, b), t);
        MM.push_back(Eigen::Triplet<float>(get_index(s.nodes, i, n),
                                           get_index(s.nodes, j, n), m));
      }
    }
  }
  return MM;
}

} // namespace femib::poisson
#endif
