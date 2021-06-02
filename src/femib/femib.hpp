#ifndef FEMIB_HPP_INCLUDED_
#define FEMIB_HPP_INCLUDED_

#include "../affine/affine.hpp"
#include "../finite_element_space/finite_element_space.hpp"
#include "../gauss/gauss.hpp"
#include "../mesh/mesh.hpp"
#include "tbb/task_scheduler_init.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <functional>
#include <tbb/concurrent_vector.h>
#include <tbb/tbb.h>
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

template <typename T> struct MeF {
  std::vector<Eigen::Triplet<T>> M;
  std::vector<Eigen::Triplet<T>> F;
};

template <typename T, int d, int e>
MeF<T>
build_diagonal(femib::finite_element_space::finite_element_space<T, d, e> s,
               femib::gauss::rule<T, d> rule,

               std::function<

                   std::function<T(femib::types::dvec<T, d>)>

                   (femib::types::F<float, 2, 1>, femib::types::F<float, 2, 1>)

                   >
                   fff) {

  tbb::task_scheduler_init init(tbb::task_scheduler_init::automatic);
  tbb::concurrent_vector<Eigen::Triplet<float>> cMM(
      s.mesh.T.size() * s.finite_element.base_functions.size() *
      s.finite_element.base_functions.size());

  tbb::concurrent_vector<Eigen::Triplet<float>> F(
      s.mesh.T.size() * s.finite_element.base_functions.size() *
      s.finite_element.base_functions.size());

  tbb::parallel_for(
      tbb::blocked_range<int>(0, s.mesh.T.size()),
      [&](const tbb::blocked_range<int> &range) {
        std::vector<Eigen::Triplet<float>> scMM;
        scMM.reserve(s.mesh.T.size() * s.finite_element.base_functions.size() *
                     s.finite_element.base_functions.size());
        std::vector<Eigen::Triplet<float>> scF;
        scF.reserve(s.mesh.T.size() * s.finite_element.base_functions.size() *
                    s.finite_element.base_functions.size());

        for (int n = range.begin(); n < range.end(); ++n) {
          femib::types::dtrian<float, 2> t = s.mesh[n];
          for (int i = 0; i < s.finite_element.base_functions.size(); ++i) {
            femib::types::F<float, 2, 1> a;

            a.dx = [&](const femib::types::dvec<float, 2> &x) {
              return (femib::affine::affineBinv(t) *
                      s.finite_element.base_functions[i].dx(
                          femib::affine::affineBinv(t) *
                          (x - femib::affine::affineb(t))));
            };

            for (int j = 0; j < s.finite_element.base_functions.size(); ++j) {

              femib::types::F<float, 2, 1> b;

              b.dx = [&](const femib::types::dvec<float, 2> &x) {
                return (femib::affine::affineBinv(t) *
                        s.finite_element.base_functions[j].dx(
                            femib::affine::affineBinv(t) *
                            (x - femib::affine::affineb(t))));
              };
              float m = femib::mesh::integrate<float, 2>(rule, fff(a, b), t);
              scMM.push_back(Eigen::Triplet<float>(
                  get_index(s.nodes, i, n), get_index(s.nodes, j, n), m));
            }

            float f_ = femib::mesh::integrate<float, 2>(
                rule,
                [&t, &s, i](femib::types::dvec<float, 2> x) {
                  float a_0 = -40 * x(0) * x(1) *
                              (s.finite_element.base_functions[i].x(
                                  femib::affine::affineBinv(t) *
                                  (x - femib::affine::affineb(t))))[0];
                  return a_0;
                },
                t);
            scF.push_back(
                Eigen::Triplet<float>(get_index(s.nodes, i, n), 0, f_));
          }
        }
        cMM.grow_by(scMM.begin(), scMM.end());
        F.grow_by(scF.begin(), scF.end());
      });

  return {{cMM.begin(), cMM.end()}, {F.begin(), F.end()}};
}

} // namespace femib::poisson
#endif
