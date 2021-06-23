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
    b1 = b1.array().min(m.P[n].array());
    b2 = b2.array().max(m.P[n].array());
  }
  box.emplace_back(b1);
  box.emplace_back(b2);
  return box;
}

template <typename T, int d>
femib::types::box<T, d> lin_spaced(const femib::types::box<T, d> &b, T delta) {
  T x_min = b[0](0);
  T x_max = b[1](0);
  T y_min = b[0](1);
  T y_max = b[1](1);
  femib::types::box<T, d> box;
  for (T x = x_min; x <= x_max; x += delta) {
    for (T y = y_min; y <= y_max; y += delta) {
      femib::types::dvec<T, d> w = {x, y};
      box.emplace_back(w);
    }
  }
  return box;
}

} // namespace femib::mesh
#endif
