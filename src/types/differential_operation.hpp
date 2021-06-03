#include "types.hpp"

template <typename T, int d, typename W>
std::function<W(femib::types::dvec<T, d>)>
operator+(const std::function<W(femib::types::dvec<T, d>)> &a,
          const std::function<W(femib::types::dvec<T, d>)> &b) {
  return ([&](const femib::types::dvec<T, d> &x) { return a(x) + b(x); });
}

template <typename T, int d, int e>
std::function<T(femib::types::dvec<T, d>)>
operator*(const std::function<T(femib::types::dvec<T, d>)> &a,
          const std::function<T(femib::types::dvec<T, d>)> &b) {
  return ([&](const femib::types::dvec<T, d> &x) { return a(x) * b(x); });
}

template <typename T, int d, int e>
std::function<T(femib::types::dvec<T, d>)>
operator*(const femib::types::F<T, d, e> &a,
          const femib::types::F<T, d, e> &b) {
  return {a.x * b.x, a.dx * b.x + b.dx * a};
}

// divergence
template <typename T, int d>
std::function<T(femib::types::dvec<T, d>)>
div(const femib::types::F<T, d, d> &A) {
  return ([&](const femib::types::dvec<T, d> &x) { return A.dx(x).trace(); });
}
