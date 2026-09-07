#ifndef P2_2d2d_HPP_INCLUDED_
#define P2_2d2d_HPP_INCLUDED_

#include "../types/types.hpp"
#include "finite_element.hpp"
#include <map>
#include <utility>
#include <vector>

namespace femib::finite_element {

template <typename T, int d, int e>
finite_element<T, d, e> create_finite_element_P2_2d2d() {
  F<T, d, e> f;
  finite_element<T, d, e> P2_2d2d;

  P2_2d2d.base_nodes = {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0},
                        {0.5, 0.5}, {0.0, 0.5}, {0.5, 0.0}};
  P2_2d2d.size = 12;

  P2_2d2d.build_nodes = [](const femib::types::mesh<T, d> &mesh) {
    femib::types::nodes<T, d> nodes;

    // assigning to each edge a sequential id and its midpoint's physical
    // location. an edge used by only one triangle is a boundary edge.
    std::map<std::pair<int, int>, int> edge_id;
    std::map<std::pair<int, int>, int> edge_count;
    std::vector<femib::types::dvec<T, d>> edge_midpoint;

    auto edge_key = [](int a, int b) {
      return a < b ? std::make_pair(a, b) : std::make_pair(b, a);
    };

    for (int n = 0; n < mesh.T.size(); ++n) {
      int v0 = mesh.T[n](0), v1 = mesh.T[n](1), v2 = mesh.T[n](2);
      std::pair<int, int> edges[3] = {edge_key(v1, v2), edge_key(v0, v2),
                                      edge_key(v0, v1)};
      for (const auto &edge : edges) {
        ++edge_count[edge];
        if (!edge_id.contains(edge)) {
          edge_id[edge] = static_cast<int>(edge_midpoint.size());
          edge_midpoint.push_back(0.5 *
                                  (mesh.P[edge.first] + mesh.P[edge.second]));
        }
      }
    }

    int size_P = mesh.P.size();
    int n_edges = static_cast<int>(edge_midpoint.size());

    nodes.P = mesh.P;
    for (femib::types::dvec<T, d> p : mesh.P) {
      nodes.P.push_back(p);
    }
    for (femib::types::dvec<T, d> m : edge_midpoint) {
      nodes.P.push_back(m);
    }
    for (femib::types::dvec<T, d> m : edge_midpoint) {
      nodes.P.push_back(m);
    }

    for (int n = 0; n < mesh.T.size(); ++n) {
      int v0 = mesh.T[n](0), v1 = mesh.T[n](1), v2 = mesh.T[n](2);
      int m12 = 2 * size_P + edge_id[edge_key(v1, v2)];
      int m02 = 2 * size_P + edge_id[edge_key(v0, v2)];
      int m01 = 2 * size_P + edge_id[edge_key(v0, v1)];

      std::vector<int> row;
      row.push_back(v0);
      row.push_back(v1);
      row.push_back(v2);
      row.push_back(m12);
      row.push_back(m02);
      row.push_back(m01);
      row.push_back(v0 + size_P);
      row.push_back(v1 + size_P);
      row.push_back(v2 + size_P);
      row.push_back(m12 + n_edges);
      row.push_back(m02 + n_edges);
      row.push_back(m01 + n_edges);
      nodes.T.push_back(row);
    }

    nodes.E = mesh.E;
    for (int i = 0; i < mesh.E.size(); ++i) {
      nodes.E.push_back(mesh.E[i] + size_P);
    }
    for (const auto &kv : edge_count) {
      if (kv.second == 1) {
        int id = edge_id[kv.first];
        nodes.E.push_back(2 * size_P + id);
        nodes.E.push_back(2 * size_P + n_edges + id);
      }
    }
    return nodes;
  };

  f.x = [](const dvec<T, d> &x) {
    if (in_std(x))
      return dvec<T, d>({1 - 3 * x(0) - 3 * x(1) + 2 * x(0) * x(0) +
                             4 * x(0) * x(1) + 2 * x(1) * x(1),
                         0});
    else
      return dvec<T, d>({0, 0});
  };
  f.dx = [](const dvec<T, d> &x) {
    if (in_std(x))
      return dmat<T, d>(
          {{-3 + 4 * x(0) + 4 * x(1), 0}, {-3 + 4 * x(0) + 4 * x(1), 0}});
    else
      return dmat<T, d>({{0, 0}, {0, 0}});
  };
  P2_2d2d.base_functions.push_back(f);

  f.x = [](const dvec<T, d> &x) {
    if (in_std(x))
      return dvec<T, d>({2 * x(0) * x(0) - x(0), 0});
    else
      return dvec<T, d>({0, 0});
  };
  f.dx = [](const dvec<T, d> &x) {
    if (in_std(x))
      return dmat<T, d>({{4 * x(0) - 1, 0}, {0, 0}});
    else
      return dmat<T, d>({{0, 0}, {0, 0}});
  };
  P2_2d2d.base_functions.push_back(f);

  f.x = [](const dvec<T, d> &x) {
    if (in_std(x))
      return dvec<T, d>({2 * x(1) * x(1) - x(1), 0});
    else
      return dvec<T, d>({0, 0});
  };
  f.dx = [](const dvec<T, d> &x) {
    if (in_std(x))
      return dmat<T, d>({{0, 0}, {4 * x(1) - 1, 0}});
    else
      return dmat<T, d>({{0, 0}, {0, 0}});
  };
  P2_2d2d.base_functions.push_back(f);

  f.x = [](const dvec<T, d> &x) {
    if (in_std(x))
      return dvec<T, d>({4 * x(0) * x(1), 0});
    else
      return dvec<T, d>({0, 0});
  };
  f.dx = [](const dvec<T, d> &x) {
    if (in_std(x))
      return dmat<T, d>({{4 * x(1), 0}, {4 * x(0), 0}});
    else
      return dmat<T, d>({{0, 0}, {0, 0}});
  };
  P2_2d2d.base_functions.push_back(f);

  f.x = [](const dvec<T, d> &x) {
    if (in_std(x))
      return dvec<T, d>({4 * x(1) - 4 * x(0) * x(1) - 4 * x(1) * x(1), 0});
    else
      return dvec<T, d>({0, 0});
  };
  f.dx = [](const dvec<T, d> &x) {
    if (in_std(x))
      return dmat<T, d>({{-4 * x(1), 0}, {4 - 4 * x(0) - 8 * x(1), 0}});
    else
      return dmat<T, d>({{0, 0}, {0, 0}});
  };
  P2_2d2d.base_functions.push_back(f);

  f.x = [](const dvec<T, d> &x) {
    if (in_std(x))
      return dvec<T, d>({4 * x(0) - 4 * x(0) * x(0) - 4 * x(0) * x(1), 0});
    else
      return dvec<T, d>({0, 0});
  };
  f.dx = [](const dvec<T, d> &x) {
    if (in_std(x))
      return dmat<T, d>({{4 - 8 * x(0) - 4 * x(1), 0}, {-4 * x(0), 0}});
    else
      return dmat<T, d>({{0, 0}, {0, 0}});
  };
  P2_2d2d.base_functions.push_back(f);

  f.x = [](const dvec<T, d> &x) {
    if (in_std(x))
      return dvec<T, d>({0, 1 - 3 * x(0) - 3 * x(1) + 2 * x(0) * x(0) +
                                4 * x(0) * x(1) + 2 * x(1) * x(1)});
    else
      return dvec<T, d>({0, 0});
  };
  f.dx = [](const dvec<T, d> &x) {
    if (in_std(x))
      return dmat<T, d>(
          {{0, -3 + 4 * x(0) + 4 * x(1)}, {0, -3 + 4 * x(0) + 4 * x(1)}});
    else
      return dmat<T, d>({{0, 0}, {0, 0}});
  };
  P2_2d2d.base_functions.push_back(f);

  f.x = [](const dvec<T, d> &x) {
    if (in_std(x))
      return dvec<T, d>({0, 2 * x(0) * x(0) - x(0)});
    else
      return dvec<T, d>({0, 0});
  };
  f.dx = [](const dvec<T, d> &x) {
    if (in_std(x))
      return dmat<T, d>({{0, 4 * x(0) - 1}, {0, 0}});
    else
      return dmat<T, d>({{0, 0}, {0, 0}});
  };
  P2_2d2d.base_functions.push_back(f);

  f.x = [](const dvec<T, d> &x) {
    if (in_std(x))
      return dvec<T, d>({0, 2 * x(1) * x(1) - x(1)});
    else
      return dvec<T, d>({0, 0});
  };
  f.dx = [](const dvec<T, d> &x) {
    if (in_std(x))
      return dmat<T, d>({{0, 0}, {0, 4 * x(1) - 1}});
    else
      return dmat<T, d>({{0, 0}, {0, 0}});
  };
  P2_2d2d.base_functions.push_back(f);

  f.x = [](const dvec<T, d> &x) {
    if (in_std(x))
      return dvec<T, d>({0, 4 * x(0) * x(1)});
    else
      return dvec<T, d>({0, 0});
  };
  f.dx = [](const dvec<T, d> &x) {
    if (in_std(x))
      return dmat<T, d>({{0, 4 * x(1)}, {0, 4 * x(0)}});
    else
      return dmat<T, d>({{0, 0}, {0, 0}});
  };
  P2_2d2d.base_functions.push_back(f);

  f.x = [](const dvec<T, d> &x) {
    if (in_std(x))
      return dvec<T, d>({0, 4 * x(1) - 4 * x(0) * x(1) - 4 * x(1) * x(1)});
    else
      return dvec<T, d>({0, 0});
  };
  f.dx = [](const dvec<T, d> &x) {
    if (in_std(x))
      return dmat<T, d>({{0, -4 * x(1)}, {0, 4 - 4 * x(0) - 8 * x(1)}});
    else
      return dmat<T, d>({{0, 0}, {0, 0}});
  };
  P2_2d2d.base_functions.push_back(f);

  f.x = [](const dvec<T, d> &x) {
    if (in_std(x))
      return dvec<T, d>({0, 4 * x(0) - 4 * x(0) * x(0) - 4 * x(0) * x(1)});
    else
      return dvec<T, d>({0, 0});
  };
  f.dx = [](const dvec<T, d> &x) {
    if (in_std(x))
      return dmat<T, d>({{0, 4 - 8 * x(0) - 4 * x(1)}, {0, -4 * x(0)}});
    else
      return dmat<T, d>({{0, 0}, {0, 0}});
  };
  P2_2d2d.base_functions.push_back(f);

  return P2_2d2d;
}
} // namespace femib::finite_element
#endif
