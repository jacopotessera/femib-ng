#ifndef P1_2d2d_HPP_INCLUDED_
#define P1_2d2d_HPP_INCLUDED_

#include "../types/types.hpp"
#include "finite_element.hpp"
#include <vector>

namespace femib::finite_element {

template <typename T, int d, int e>
finite_element<T, d, e> create_finite_element_P1_2d2d() {
  F<T, d, e> f;
  finite_element<T, d, e> P1_2d2d;

  P1_2d2d.base_nodes = {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}};
  P1_2d2d.size = 6;
  P1_2d2d.build_nodes = [](const femib::types::mesh<T, d> &mesh) {
    femib::types::nodes<T, d> nodes;
    nodes.P = mesh.P;
    for (int n = 0; n < mesh.T.size(); ++n) {
      std::vector<int> row;
      nodes.T.push_back(row);
      for (int i = 0; i < mesh.T[n].size() * d; ++i) {
        nodes.T[n].push_back(mesh.T[n](i % d));
      }
    }
    nodes.E = mesh.E;
    return nodes;
  };

  f.x = [](const dvec<T, d> &x) {
    if (in_std(x))
      return dvec<T, d>({1 - x(0) - x(1), 0});
    else
      return dvec<T, d>({0, 0});
  };

  f.dx = [](const dvec<T, d> &x) {
    if (in_std(x))
      return dmat<T, d>({{-1, -1}, {0, 0}});
    else
      return dmat<T, d>({{0, 0}, {0, 0}});
  };

  P1_2d2d.base_functions.push_back(f);

  f.x = [](const dvec<T, d> &x) {
    if (in_std(x))
      return dvec<T, d>({x(0), 0});
    else
      return dvec<T, d>({0, 0});
  };

  f.dx = [](const dvec<T, d> &x) {
    if (in_std(x))
      return dmat<T, d>({{1, 0}, {0, 0}});
    else
      return dmat<T, d>({{0, 0}, {0, 0}});
  };

  P1_2d2d.base_functions.push_back(f);

  f.x = [](const dvec<T, d> &x) {
    if (in_std(x))
      return dvec<T, d>({x(1), 0});
    else
      return dvec<T, d>({0, 0});
  };

  f.dx = [](const dvec<T, d> &x) {
    if (in_std(x))
      return dmat<T, d>({{0, 1}, {0, 0}});
    else
      return dmat<T, d>({{0, 0}, {0, 0}});
  };

  P1_2d2d.base_functions.push_back(f);

  f.x = [](const dvec<T, d> &x) {
    if (in_std(x))
      return dvec<T, d>({0, 1 - x(0) - x(1)});
    else
      return dvec<T, d>({0, 0});
  };

  f.dx = [](const dvec<T, d> &x) {
    if (in_std(x))
      return dmat<T, d>({{0, 0}, {-1, -1}});
    else
      return dmat<T, d>({{0, 0}, {0, 0}});
  };

  P1_2d2d.base_functions.push_back(f);

  f.x = [](const dvec<T, d> &x) {
    if (in_std(x))
      return dvec<T, d>({0, x(0)});
    else
      return dvec<T, d>({0, 0});
  };

  f.dx = [](const dvec<T, d> &x) {
    if (in_std(x))
      return dmat<T, d>({{0, 0}, {1, 0}});
    else
      return dmat<T, d>({{0, 0}, {0, 0}});
  };

  P1_2d2d.base_functions.push_back(f);

  f.x = [](const dvec<T, d> &x) {
    if (in_std(x))
      return dvec<T, d>({0, x(1)});
    else
      return dvec<T, d>({0, 0});
  };

  f.dx = [](const dvec<T, d> &x) {
    if (in_std(x))
      return dmat<T, d>({{0, 0}, {0, 1}});
    else
      return dmat<T, d>({{0, 0}, {0, 0}});
  };

  P1_2d2d.base_functions.push_back(f);

  return P1_2d2d;
}
} // namespace femib::finite_element
#endif
