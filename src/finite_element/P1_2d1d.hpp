#ifndef P1_2d1d_HPP_INCLUDED_
#define P1_2d1d_HPP_INCLUDED_

#include "../types/types.hpp"
#include "finite_element.hpp"
#include <vector>

namespace femib::finite_element {

template <typename T, int d, int e>
finite_element<T, d, e> create_finite_element_P1_2d1d() {
  function_with_gradient<T, d, e> f;
  finite_element<T, d, e> P1_2d1d;

  P1_2d1d.base_nodes = {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}};
  P1_2d1d.size = 3;
  P1_2d1d.build_nodes = [](const femib::types::mesh<T, d> &mesh) {
    femib::types::nodes<T, d> nodes;
    nodes.P = mesh.P;
    for (int n = 0; n < std::ssize(mesh.T); ++n) {
      std::vector<int> row;
      nodes.T.push_back(row);
      for (int i = 0; i < mesh.T[n].size(); ++i) {
        nodes.T[n].push_back(mesh.T[n](i));
      }
    }
    nodes.E = mesh.E;
    return nodes;
  };

  f.value = [](const dvec<T, d> &x) {
    if (is_point_in_reference_triangle(x))
      return (dvec<T, e>){1 - x(0) - x(1)};
    else
      return (dvec<T, e>){0.0};
  };

  f.gradient = [](const dvec<T, d> &x) {
    if (is_point_in_reference_triangle(x))
      return (rmat<T, d, e>){{-1, -1}};
    else
      return (rmat<T, d, e>){{0, 0}};
  };

  P1_2d1d.base_functions.push_back(f);

  f.value = [](const dvec<T, d> &x) {
    if (is_point_in_reference_triangle(x))
      return (dvec<T, e>){x(0)};
    else
      return (dvec<T, e>){0.0};
  };

  f.gradient = [](const dvec<T, d> &x) {
    if (is_point_in_reference_triangle(x))
      return (rmat<T, d, e>){{1, 0}};
    else
      return (rmat<T, d, e>){{0, 0}};
  };

  P1_2d1d.base_functions.push_back(f);

  f.value = [](const dvec<T, d> &x) {
    if (is_point_in_reference_triangle(x))
      return (dvec<T, e>){x(1)};
    else
      return (dvec<T, e>){0.0};
  };

  f.gradient = [](const dvec<T, d> &x) {
    if (is_point_in_reference_triangle(x))
      return (rmat<T, d, e>){{0, 1}};
    else
      return (rmat<T, d, e>){{0, 0}};
  };

  P1_2d1d.base_functions.push_back(f);

  return P1_2d1d;
}

} // namespace femib::finite_element
#endif
