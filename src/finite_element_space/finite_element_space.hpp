#ifndef FINITE_ELEMENT_SPACE_HPP_INCLUDED_
#define FINITE_ELEMENT_SPACE_HPP_INCLUDED_
#include "../affine/affine.hpp"
#include "../cuda/cuda.h"
#include "../finite_element/finite_element.hpp"
#include "../mesh/mesh.hpp"
#include "../types/types.hpp"
#include <ranges>

namespace femib::finite_element_space {

template <typename T, int d, int e> struct finite_element_space {
  femib::finite_element::finite_element<T, d, e> finite_element;
  femib::types::mesh<T, d> mesh;
  femib::types::nodes<T, d> nodes;

  std::vector<femib::types::dvec<T, e>>
  interpolate(const Eigen::Matrix<T, Eigen::Dynamic, 1> &u,
              std::vector<femib::types::dvec<T, d>> &points) {
    std::vector<int> points_positions =
        femib::mesh::find_points<T, d>(mesh, points);

    std::vector<femib::types::dvec<T, e>> result(
        points.size(), femib::types::dvec<T, e>::Zero());
    for (size_t k = 0; k < points.size(); ++k) {
      if (points_positions[k] < 0) {
        continue;
      }
      femib::types::dtrian<T, d> t = mesh[points_positions[k]];
      femib::types::dvec<T, d> ref =
          femib::affine::affine_inv<T, d>(t, points[k]);
      for (int j = 0; j < finite_element.base_functions.size(); ++j) {
        femib::types::F<T, d, e> phi = finite_element.base_functions[j];
        result[k] += u(nodes.get_index(j, points_positions[k])) * phi.x(ref);
      }
    }
    return result;
  }

  std::vector<std::pair<types::dvec<T, d>, types::dvec<T, e>>>
  plot(Eigen::Matrix<T, Eigen::Dynamic, 1> u, T delta) {
    femib::types::box<T, d> box =
        femib::mesh::lin_spaced<T, d>(femib::mesh::find_box<T, d>(mesh), delta);
    std::vector<types::dvec<T, e>> results = interpolate(u, box);

    std::vector<std::pair<types::dvec<T, d>, types::dvec<T, e>>> plot_data;
    for (auto [p, r] : std::views::zip(box, results)) {
      plot_data.emplace_back(
          std::pair<types::dvec<T, d>, types::dvec<T, e>>(p, r));
    }
    return plot_data;
  }
};
} // namespace femib::finite_element_space
#endif
