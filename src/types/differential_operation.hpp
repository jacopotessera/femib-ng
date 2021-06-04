#include "types.hpp"

template <typename T, int d>
std::function<T(femib::types::dvec<T, d>)>
operator*(const std::function<T(femib::types::dvec<T, d>)> &a,
          const std::function<T(femib::types::dvec<T, d>)> &b) {
  return ([&](const femib::types::dvec<T, d> &x) { return a(x) * b(x); });
}

// divergence
template <typename T, int d>
std::function<T(femib::types::dvec<T, d>)>
div(const femib::types::F<T, d, d> &A) {
  return [&](const femib::types::dvec<T, d> &x) { return A.dx(x).trace(); };
}

template <typename T, int d>
std::function<T(femib::types::dvec<T, d>)> project(
    const std::function<femib::types::dvec<T, d>(femib::types::dvec<T, d>)> &a,
    int i) {
  return ([&, i](const femib::types::dvec<T, d> &x) { return a(x)(i); });
}
