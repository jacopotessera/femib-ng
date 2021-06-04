#ifndef FINITE_ELEMENT_HPP_INCLUDED_
#define FINITE_ELEMENT_HPP_INCLUDED_

#include "../types/types.hpp"
#include <vector>

namespace femib::finite_element {

template <typename T, int d> using dvec = femib::types::dvec<T, d>;
template <typename T, int d> using dmat = femib::types::dmat<T, d>;
template <typename T, int d, int e> using rmat = femib::types::rmat<T, d, e>;
template <typename T, int d, int e> using F = femib::types::F<T, d, e>;
template <typename T, int d> using mesh_ = femib::types::mesh<T, d>;

template <typename T, int d, int e> struct finite_element {
  std::vector<F<T, d, e>> base_functions;
  std::vector<dvec<T, d>> base_nodes;
  femib::types::nodes<T, d> nodes;
  std::function<femib::types::nodes<T, d>(const femib::types::mesh<T, d>)>
      build_nodes;
  int size;

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

} // namespace femib::finite_element
#endif
