#ifndef MESH_HPP_INCLUDED_
#define MESH_HPP_INCLUDED_

#include "../gauss/gauss.hpp"
#include "../types/types.hpp"
#include "../read/read.hpp"
#include "../affine/affine.hpp"
#include <string>

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
T femib::mesh::integrate(const femib::gauss::rule<T, d> &rule,
                         const std::function<T(femib::types::dvec<T, d>)> &f,
                         const femib::types::mesh<T, d> &mesh) {
  T integral = 0.0;
  for (const femib::types::dtrian<T, d> &t : mesh) {
    integral += femib::mesh::integrate(rule, f, t);
  }
  return integral;
}

template <typename T, int d>
femib::types::mesh<T, d> femib::mesh::read(const std::string &filename_p,
                                           const std::string &filename_t,
                                           const std::string &filename_e) {
  femib::types::mesh<T, d> mesh = read_mesh_file<T,d>(filename_p, filename_t, filename_e);
  return mesh;
}

} // namespace femib::mesh
#endif
