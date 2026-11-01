#ifndef FEMIB_AFFINE_HPP
#define FEMIB_AFFINE_HPP

#include <Eigen/Core>
#include <Eigen/Dense>
#include <cmath>

#include "../types/types.hpp"

// The affine map x = B x_ref + b from the reference triangle (0,0), (1,0),
// (0,1) onto the triangle t: b is its translation, B its jacobian.
namespace femib::affine {

template <typename T, int d>
femib::types::dvec<T, d> translation(const femib::types::dtrian<T, d> &t) {
  return t[0];
}

template <typename T, int d>
femib::types::dmat<T, d> jacobian(const femib::types::dtrian<T, d> &t) {
  femib::types::dmat<T, d> B;
  if (d == 1)
    B << t[1](0) - t[0](0);
  if (d == 2)
    B << t[1](0) - t[0](0), t[2](0) - t[0](0), t[1](1) - t[0](1),
        t[2](1) - t[0](1);
  if (d == 3)
    B << t[1](0) - t[0](0), t[2](0) - t[0](0), t[3](0) - t[0](0),
        t[1](1) - t[0](1), t[2](1) - t[0](1), t[3](1) - t[0](1),
        t[1](2) - t[0](2), t[2](2) - t[0](2), t[3](2) - t[0](2);
  return B;
}

template <typename T, int d>
femib::types::dmat<T, d> inverse_jacobian(const femib::types::dtrian<T, d> &t) {
  return jacobian(t).inverse();
}

template <typename T, int d>
T abs_jacobian_determinant(const femib::types::dtrian<T, d> &t) {
  return std::abs(jacobian(t).determinant());
}

template <typename T, int d>
femib::types::dvec<T, d>
reference_to_physical(const femib::types::dtrian<T, d> &t,
                      const femib::types::dvec<T, d> &x) {
  return jacobian(t) * x + translation(t);
}

template <typename T, int d>
femib::types::dvec<T, d>
physical_to_reference(const femib::types::dtrian<T, d> &t,
                      const femib::types::dvec<T, d> &x) {
  return jacobian(t).inverse() * (x - translation(t));
}

} // namespace femib::affine
#endif
