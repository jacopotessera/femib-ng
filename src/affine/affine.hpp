#ifndef FEMIB_AFFINE_HPP
#define FEMIB_AFFINE_HPP

#include <Eigen/Core>
#include <Eigen/Dense>
#include <cmath>

#include "../types/types.hpp"

namespace femib::affine {

template <typename T, int d>
femib::types::dvec<T, d> affineb(const femib::types::dtrian<T, d> &t) {
  return t[0];
}

template <typename T, int d>
femib::types::dmat<T, d> affineB(const femib::types::dtrian<T, d> &t) {
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
femib::types::dmat<T, d> affineBinv(const femib::types::dtrian<T, d> &t) {
  return affineB(t).inverse();
}

template <typename T, int d> T affineBdet(const femib::types::dtrian<T, d> &t) {
  return abs(affineB(t).determinant());
}

template <typename T, int d>
femib::types::dvec<T, d> affine(const femib::types::dtrian<T, d> &t,
                                const femib::types::dvec<T, d> &x) {
  return affineB(t) * x + affineb(t);
}

template <typename T, int d>
femib::types::dvec<T, d> affine_inv(const femib::types::dtrian<T, d> &t,
                                    const femib::types::dvec<T, d> &x) {
  return affineB(t).inverse() * (x - affineb(t));
}

} // namespace femib::affine
#endif
