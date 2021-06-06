#ifndef P1_B_2d2d_HPP_INCLUDED_
#define P1_B_2d2d_HPP_INCLUDED_

#include "../types/types.hpp"
#include "finite_element.hpp"
#include <vector>

namespace femib::finite_element {

template <typename T, int d>
dvec<T, d> find_center_of(const femib::types::dtrian<T, d> &t) {
  dvec<T, d> c = {0.0, 0.0};
  double k = 1.0 / t.size();
  for (int j = 0; j < t.size(); ++j) {
    c = c + k * t[j];
  }
  return c;
}

template <typename T, int d, int e>
finite_element<T, d, e> create_finite_element_P1_B_2d2d() {
  F<T, d, e> f;
  finite_element<T, d, e> P1_2d2d;

  P1_2d2d.base_nodes = {
      {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}, {1.0 / 3.0, 1.0 / 3.0}};
  P1_2d2d.size = 8;
  P1_2d2d.build_nodes = [](const femib::types::mesh<T, d> &mesh) {
    femib::types::nodes<T, d> nodes;
    nodes.P = mesh.P;
    for (femib::types::dvec<T, d> p : mesh.P) {
      nodes.P.push_back(p);
    }

    for (int n = 0; n < mesh.T.size(); ++n) {
      femib::types::dtrian<T, d> t = mesh[n];
      nodes.P.push_back(find_center_of(t));
    }
    for (int n = 0; n < mesh.T.size(); ++n) {
      femib::types::dtrian<T, d> t = mesh[n];
      nodes.P.push_back(find_center_of(t));
    }

    int size_P = mesh.P.size();
    int size_T = mesh.T.size();
    for (int n = 0; n < mesh.T.size(); ++n) {
      std::vector<int> row;
      nodes.T.push_back(row);
      for (int i = 0; i < mesh.T[n].size(); ++i) {
        nodes.T[n].push_back(mesh.T[n](i));
      }
      for (int i = 0; i < mesh.T[n].size(); ++i) {
        nodes.T[n].push_back(mesh.T[n](i) + size_P);
      }
      nodes.T[n].push_back(n + 2 * size_P);
      nodes.T[n].push_back(n + size_T + 2 * size_P);
    }

    nodes.E = mesh.E;
    for (int i = 0; i < mesh.E.size(); ++i) {
      nodes.E.push_back(mesh.E[i] + size_P);
    }
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

  f.x = [](const dvec<T, d> &x) {
    if (in_std(x))
      return dvec<T, d>({x(0) * x(0) + 7 * x(0) * x(1) + x(1) * x(1), 0});
    else
      return dvec<T, d>({0, 0});
  };

  f.dx = [](const dvec<T, d> &x) {
    if (in_std(x))
      return dmat<T, d>({{2 * x(0) + 7 * x(1), 2 * x(1) + 7 * x(0)}, {0, 0}});
    else
      return dmat<T, d>({{0, 0}, {0, 0}});
  };

  P1_2d2d.base_functions.push_back(f);

  f.x = [](const dvec<T, d> &x) {
    if (in_std(x))
      return dvec<T, d>({0, x(0) * x(0) + 7 * x(0) * x(1) + x(1) * x(1)});
    else
      return dvec<T, d>({0, 0});
  };

  f.dx = [](const dvec<T, d> &x) {
    if (in_std(x))
      return dmat<T, d>({{0, 0}, {2 * x(0) + 7 * x(1), 2 * x(1) + 7 * x(0)}});
    else
      return dmat<T, d>({{0, 0}, {0, 0}});
  };

  P1_2d2d.base_functions.push_back(f);

  return P1_2d2d;
}
} // namespace femib::finite_element
#endif
