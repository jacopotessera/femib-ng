#ifndef OPERATORS_HPP_INCLUDED_
#define OPERATORS_HPP_INCLUDED_
#include "types.hpp"

namespace femib::operators {

template <typename T, int d>
std::function<T(femib::types::dvec<T, d>)>
divergence(const femib::types::function_with_gradient<T, d, d> &a) {
  return
      [=](const femib::types::dvec<T, d> &x) { return a.gradient(x).trace(); };
}

template <typename T, int d>
std::function<femib::types::dmat<T, d>(femib::types::dvec<T, d>)>
symmetric_gradient(const femib::types::function_with_gradient<T, d, d> &a) {
  return ([=](const femib::types::dvec<T, d> &x) -> femib::types::dmat<T, d> {
    return T(0.5) * (a.gradient(x) + a.gradient(x).transpose());
  });
}

template <typename T, int d>
T double_dot(const femib::types::dmat<T, d> &A,
             const femib::types::dmat<T, d> &B) {
  return (A.transpose() * B).trace();
}

template <typename T, int d>
std::function<T(femib::types::dvec<T, d>)> symmetric_gradient_double_dot(
    const femib::types::function_with_gradient<T, d, d> &a,
    const femib::types::function_with_gradient<T, d, d> &b) {
  return ([=](const femib::types::dvec<T, d> &x) {
    return double_dot(symmetric_gradient(a)(x), symmetric_gradient(b)(x));
  });
}

template <typename T, int d>
std::function<femib::types::dvec<T, d>(femib::types::dvec<T, d>)>
convection(const femib::types::function_with_gradient<T, d, d> &w,
           const femib::types::function_with_gradient<T, d, d> &u) {
  return [=](const femib::types::dvec<T, d> &x) -> femib::types::dvec<T, d> {
    // ((w . grad) u) (x) =  u.gradient(x).transpose() * w.value(x)
    return u.gradient(x).transpose() * w.value(x);
  };
}

template <typename T, int d, int e>
std::function<T(femib::types::dvec<T, d>)>
mass(femib::types::function_with_gradient<T, d, e> a,
     femib::types::function_with_gradient<T, d, e> b) {
  return [a, b](const femib::types::dvec<T, d> &x) {
    return a.value(x).dot(b.value(x));
  };
}

template <typename T, int d, int e>
std::function<T(femib::types::dvec<T, d>)>
gradient_double_dot(femib::types::function_with_gradient<T, d, e> a,
                    femib::types::function_with_gradient<T, d, e> b) {
  return [a, b](const femib::types::dvec<T, d> &x) {
    return a.gradient(x).cwiseProduct(b.gradient(x)).sum();
  };
}

} // namespace femib::operators
#endif
