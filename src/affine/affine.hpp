#ifndef FEMIB_AFFINE_HPP
#define FEMIB_AFFINE_HPP

#include <functional>

#include "../types/types.hpp"
#include <Eigen/Core>

template <typename T, int d>
femib::types::dvec<T, d> affine(const femib::types::dtrian<T, d> &t,
                                const femib::types::dvec<T, d> &x);
template <typename T, int d>
femib::types::dvec<T, d> affine_inv(const femib::types::dtrian<T, d> &t,
                                    const femib::types::dvec<T, d> &x);
template <typename T, int d>
femib::types::dmat<T, d> affineB(const femib::types::dtrian<T, d> &t);
template <typename T, int d>
femib::types::dvec<T, d> affineb(const femib::types::dtrian<T, d> &t);

#endif
