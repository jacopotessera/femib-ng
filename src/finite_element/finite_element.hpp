#ifndef FINITE_ELEMENT_HPP_INCLUDED_
#define FINITE_ELEMENT_HPP_INCLUDED_

#include "../types/types.hpp"
#include <vector>

namespace femib::finite_element {

template <typename T, int d> using dvec = femib::types::dvec<T, d>;
template <typename T, int d> using dmat = femib::types::dmat<T, d>;
template <typename T, int d> using F = femib::types::F<T, d>;

template <typename T, int d, int e> struct finite_element {
  std::vector<F<T, d>> base_functions;
  std::vector<dvec<T, d>> base_nodes;
  int dim1 = d;
  int dim2 = e;
};

template <typename T, int d> bool in_std(const dvec<T, d> &P) {
  bool e = true;
  T q = 0;
  for (T p : P) {
    e = e && p >= 0;
    q += p;
  }
  return e && q <= 1;
};

template <typename T, int d, int e>
finite_element<T, d, e> create_finite_element_P1_2d1d() {
  F<T, d> f;
  finite_element<T, d, e> P1_2d1d;

  P1_2d1d.base_nodes = {{0.0,0.0},{1.0,0.0},{0.0,1.0}};
  f.x = [](const dvec<T, d> &x) {
    if (in_std(x))
      return (dvec<T,d>){1 - x(0) - x(1), 0.0};
    else
      return (dvec<T,d>){0.0, 0.0};
  };

  f.dx = [](const dvec<T, d> &x) {
    if (in_std(x))
      return (dmat<T, d>){{-1, -1}, {0, 0}};
    else
      return (dmat<T, d>){{0, 0}, {0, 0}};
  };

  P1_2d1d.base_functions.push_back(f);

  f.x = [](const dvec<T, d> &x) {
    if (in_std(x))
      return (dvec<T, d>){x(0), 0.0};
    else
      return (dvec<T, d>){0.0, 0.0};
  };

  f.dx = [](const dvec<T, d> &x) {
    if (in_std(x))
      return (dmat<T, d>){{1, 0}, {0, 0}};
    else
      return (dmat<T, d>){{0, 0}, {0, 0}};
  };

  P1_2d1d.base_functions.push_back(f);

  f.x = [](const dvec<T, d> &x) {
    if (in_std(x))
      return (dvec<T, d>){x(1), 0.0};
    else
      return (dvec<T, d>){0.0, 0.0};
  };

  f.dx = [](const dvec<T, d> &x) {
    if (in_std(x))
      return (dmat<T, d>){{0, 1}, {0, 0}};
    else
      return (dmat<T, d>){{0, 0}, {0, 0}};
  };

  P1_2d1d.base_functions.push_back(f);

  return P1_2d1d;
}

} // namespace femib::finite_element
#endif
