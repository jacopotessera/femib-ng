#include "mesh.hpp"
#include "../affine/affine.hpp"
#include "spdlog/spdlog.h"

template <typename T, int d>
T femib::mesh::integrate(const femib::gauss::rule<T, d> &rule,
                         const std::function<T(femib::types::dvec<T, d>)> &f,
                         const femib::types::dtrian<T, d> &t) {

  std::function<T(femib::types::dvec<T, d>)> g =
      [&, f](const femib::types::dvec<T, d> &x) {
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

template float femib::mesh::integrate<float,2>(
    const femib::gauss::rule<float, 2> &rule,
    const std::function<float(femib::types::dvec<float, 2>)> &f,
    const femib::types::mesh<float, 2> &mesh);
