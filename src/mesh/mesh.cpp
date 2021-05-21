#include "mesh.hpp"
#include "../affine/affine.hpp"
#include "spdlog/spdlog.h"

#include "read.cpp"

template <typename T, int d>
T femib::mesh::integrate(const femib::gauss::rule<T, d> &rule,
                         const std::function<T(femib::types::dvec<T, d>)> &f,
                         const femib::types::dtrian<T, d> &t) {

  std::function<T(femib::types::dvec<T, d>)> g =
      [&t, &f](const femib::types::dvec<T, d> &x) {
        spdlog::debug("[mesh] node x  {}", x[0]);
        spdlog::debug("[mesh] node y  {}", x[1]);
        spdlog::debug("[mesh] f(node) {}", f(affine(t, x)));
        spdlog::debug("[mesh] f(0) {}", f({0.0, 0.0}));
        return affineBdet(t) * f(affine(t, x));
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
  femib::types::mesh<T, d> mesh = readMesh(filename_p, filename_t, filename_e);
  return mesh;
}

template float femib::mesh::integrate<float, 2>(
    const femib::gauss::rule<float, 2> &rule,
    const std::function<float(femib::types::dvec<float, 2>)> &f,
    const femib::types::mesh<float, 2> &mesh);
template femib::types::mesh<float, 2>
femib::mesh::read<float, 2>(const std::string &filename_p,
                            const std::string &filename_t,
                            const std::string &filename_e);
