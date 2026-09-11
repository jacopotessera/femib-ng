#ifndef DIFFERENTIAL_OPERATION_HPP_INCLUDED_
#define DIFFERENTIAL_OPERATION_HPP_INCLUDED_

#include "types.hpp"

template <typename T, int d>
std::function<T(femib::types::dvec<T, d>)>
operator*(const std::function<T(femib::types::dvec<T, d>)> &a,
          const std::function<T(femib::types::dvec<T, d>)> &b) {
  return ([=](const femib::types::dvec<T, d> &x) { return a(x) * b(x); });
}

// divergence
template <typename T, int d>
std::function<T(femib::types::dvec<T, d>)>
div(const femib::types::F<T, d, d> &a) {
  return [=](const femib::types::dvec<T, d> &x) { return a.dx(x).trace(); };
}

template <typename T, int d, int e>
std::function<T(femib::types::dvec<T, d>)> project(
    const std::function<femib::types::dvec<T, e>(femib::types::dvec<T, d>)> &a,
    int i) {
  return ([=](const femib::types::dvec<T, d> &x) { return a(x)(i); });
}

template <typename T, int d>
std::function<femib::types::dmat<T, d>(femib::types::dvec<T, d>)>
symm(const femib::types::F<T, d, d> &a) {
  return ([=](const femib::types::dvec<T, d> &x) -> femib::types::dmat<T, d> {
    return T(0.5) * (a.dx(x) + a.dx(x).transpose());
  });
}

// doppio prodotto interno
template <typename T, int d>
T dpi(const femib::types::dmat<T, d> &A, const femib::types::dmat<T, d> &B) {
  return (A.transpose() * B).trace();
}

template <typename T, int d>
std::function<T(femib::types::dvec<T, d>)>
dpi(const femib::types::F<T, d, d> &a, const femib::types::F<T, d, d> &b) {
  return ([=](const femib::types::dvec<T, d> &x) {
    return dpi(symm(a)(x), symm(b)(x));
  });
}

template <typename T, int d>
std::function<femib::types::dvec<T, d>(femib::types::dvec<T, d>)>
convection(const femib::types::F<T, d, d> &w,
           const femib::types::F<T, d, d> &u) {
  return [=](const femib::types::dvec<T, d> &x) -> femib::types::dvec<T, d> {
    // ((w . grad) u) (x) =  u.dx(x).transpose() * w.x(x)
    return u.dx(x).transpose() * w.x(x);
  };
}

template <typename T, int d>
std::function<T(femib::types::dvec<T, d>)> mass(femib::types::F<T, d, d> a,
                                                femib::types::F<T, d, d> b) {
  return
      [a, b](const femib::types::dvec<T, d> &x) { return a.x(x).dot(b.x(x)); };
}

template <typename T, int d>
std::function<T(femib::types::dvec<T, d>)> zero(femib::types::F<T, d, d>) {
  return [](const femib::types::dvec<T, d> &) { return T(0); };
}

#endif
