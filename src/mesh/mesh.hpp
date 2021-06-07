#ifndef MESH_HPP_INCLUDED_
#define MESH_HPP_INCLUDED_

#include "../affine/affine.hpp"
#include "../gauss/gauss.hpp"
#include "../read/read.hpp"
#include "../types/types.hpp"
#include <execution>
#include <functional>
#include <numeric>
#include <string>

namespace femib::mesh {

template <typename T, int d>
T integrate(const femib::gauss::rule<T, d> &rule,
            const std::function<T(femib::types::dvec<T, d>)> &f,
            const femib::types::dtrian<T, d> &t) {

  std::function<T(femib::types::dvec<T, d>)> g =
      [&t, &f](const femib::types::dvec<T, d> &x) {
        return femib::affine::affineBdet(t) * f(femib::affine::affine(t, x));
      };
  return femib::gauss::integrate<T, d>(rule, g);
}

template <typename T, int d>
T integrate(const femib::gauss::rule<T, d> &rule,
            const std::function<T(femib::types::dvec<T, d>)> &f,
            const femib::types::mesh<T, d> &mesh) {
  auto unary_op = [&rule, &f](const femib::types::dtrian<T, d> &t) {
    return femib::mesh::integrate(rule, f, t);
  };
  return std::transform_reduce(std::execution::seq, mesh.N.begin(),
                               mesh.N.end(), 0.0, std::plus<>(), unary_op);
}

template <typename T, int d>
femib::types::mesh<T, d> read(const std::string &filename_p,
                              const std::string &filename_t,
                              const std::string &filename_e) {
  femib::types::mesh<T, d> mesh =
      read_mesh_file<T, d>(filename_p, filename_t, filename_e);
  return mesh;
}

template <typename T, int d>
femib::types::box<T, d> find_box(const femib::types::mesh<T, d> &m) {
  femib::types::box<T, d> box;
  femib::types::dvec<T, d> b1 = m.P[0];
  femib::types::dvec<T, d> b2 = m.P[0];
  for (int n = 1; n < m.P.size(); ++n) {
    femib::types::dvec<T, d> v = m.P[n];

    if (b1(0) > v(0))
      b1(0) = v(0);
    if (b1(1) > v(1))
      b1(1) = v(1);

    if (b2(0) < v(0))
      b2(0) = v(0);
    if (b2(1) < v(1))
      b2(1) = v(1);
  }
  box.emplace_back(b1);
  box.emplace_back(b2);
  return box;
}

} // namespace femib::mesh
#endif
