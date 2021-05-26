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
#include <vector>

namespace femib::mesh {

template <typename T, int d>
T integrate(const femib::gauss::rule<T, d> &rule,
            const std::function<T(femib::types::dvec<T, d>)> &f,
            const femib::types::dtrian<T, d> &t);
template <typename T, int d>
T integrate(const femib::gauss::rule<T, d> &rule,
            const std::function<T(femib::types::dvec<T, d>)> &f,
            const femib::types::mesh<T, d> &mesh);
template <typename T, int d>
femib::types::mesh<T, d> read(const std::string &filename_p,
                              const std::string &filename_t,
                              const std::string &filename_e);

template <typename T, int d>
std::vector<femib::types::dtrian<T, d>>
init(const femib::types::mesh<T, d> &mesh);

template <typename T, int d>
T femib::mesh::integrate(const femib::gauss::rule<T, d> &rule,
                         const std::function<T(femib::types::dvec<T, d>)> &f,
                         const femib::types::dtrian<T, d> &t) {

  std::function<T(femib::types::dvec<T, d>)> g =
      [&t, &f](const femib::types::dvec<T, d> &x) {
        return femib::affine::affineBdet(t) * f(femib::affine::affine(t, x));
      };
  return femib::gauss::integrate<T, d>(rule, g);
}

template <typename T, int d>
std::vector<femib::types::dtrian<T, d>>
femib::mesh::init(const femib::types::mesh<T, d> &mesh) {
  std::vector<femib::types::dtrian<T, d>> dtrianfd;
  for (femib::types::ditrian<d> triangle : mesh.T) {
    femib::types::dtrian<T, d> t;
    t.reserve(d + 1);
    for (int point : triangle) {
      t.emplace_back(mesh.P[point]);
    }
    dtrianfd.emplace_back(t);
  }
  return dtrianfd;
}

template <typename T, int d>
T femib::mesh::integrate(const femib::gauss::rule<T, d> &rule,
                         const std::function<T(femib::types::dvec<T, d>)> &f,
                         const femib::types::mesh<T, d> &mesh) {
  auto unary_op = [&rule, &f](const femib::types::dtrian<T, d> &t) {
    return femib::mesh::integrate(rule, f, t);
  };
  std::vector<femib::types::dtrian<T, d>> m = femib::mesh::init(mesh);
  return std::transform_reduce(std::execution::seq, m.begin(), m.end(), 0.0,
                               std::plus<>(), unary_op);
}

template <typename T, int d>
femib::types::mesh<T, d> femib::mesh::read(const std::string &filename_p,
                                           const std::string &filename_t,
                                           const std::string &filename_e) {
  femib::types::mesh<T, d> mesh =
      read_mesh_file<T, d>(filename_p, filename_t, filename_e);
  return mesh;
}

} // namespace femib::mesh
#endif
