#ifndef P0_2d1d_HPP_INCLUDED_
#define P0_2d1d_HPP_INCLUDED_

#include "../types/types.hpp"
#include "finite_element.hpp"
#include <vector>

namespace femib::finite_element {

template <typename T, int d>
dvec<T, d> find_center_of(const femib::types::mesh<T, d> &mesh, int n) {
  dvec<T, d> c = {0.0, 0.0};
  double k = 1.0 / mesh.T[n].size();
  for (int j = 0; j < mesh.T[n].size(); ++j) {
    c = c + k * mesh.P[mesh.T[n](j)];
  }
  return c;
}

template <typename T, int d, int e>
finite_element<T, d, e> create_finite_element_P0_2d1d() {
  F<T, d, e> f;
  finite_element<T, d, e> P0_2d1d;

  P0_2d1d.base_nodes = {{1.0 / 3.0, 1.0 / 3.0}};
  P0_2d1d.size = 1;
  P0_2d1d.build_nodes = [](const femib::types::mesh<T, d> &mesh) {
    femib::types::nodes<T, d> nodes;
    for (int n = 0; n < mesh.T.size(); ++n) {
      std::vector<int> row;
      nodes.T.push_back(row);
      nodes.P.push_back(find_center_of<T, d>(mesh, n));
      nodes.T[n].push_back(n);
    }
    nodes.E = mesh.E;
    return nodes;
  };

  f.x = [](const dvec<T, d> &x) {
    if (in_std(x))
      return (dvec<T, e>){1.0};
    else
      return (dvec<T, e>){0.0};
  };

  f.dx = [](const dvec<T, d> &x) { return (rmat<T, d, e>){{0, 0}}; };

  P0_2d1d.base_functions.push_back(f);

  return P0_2d1d;
}

} // namespace femib::finite_element
#endif
