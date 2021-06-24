#ifndef FINITE_ELEMENT_SPACE_HPP_INCLUDED_
#define FINITE_ELEMENT_SPACE_HPP_INCLUDED_
#include "../affine/affine.hpp"
#include "../cuda/cuda.h"
#include "../finite_element/finite_element.hpp"
#include "../mesh/mesh.hpp"
#include "../types/types.hpp"

namespace femib::finite_element_space {

template <typename T, int d, int e> struct finite_element_space {
  femib::finite_element::finite_element<T, d, e> finite_element;
  femib::types::mesh<T, d> mesh;
  femib::types::nodes<T, d> nodes;

  std::vector<std::vector<std::vector<T>>>
  plot(Eigen::Matrix<T, Eigen::Dynamic, 1> xx) {
    std::vector<std::vector<std::vector<float>>> uuu;

    femib::types::box<T, d> box = femib::mesh::find_box<T, d>(mesh);

    femib::types::box<T, d> boxx = femib::mesh::lin_spaced<T, d>(box, 0.027);
    bool N[boxx.size() * mesh.N.size()];

    femib::cuda::serial_accurate<T, d>(boxx.data(), boxx.size(), mesh.N.data(),
                                       mesh.N.size(), N);

    std::vector<int> NNN;

    for (int i = 0; i < boxx.size(); ++i) {
      for (int n = 0; n < mesh.N.size(); ++n) {
        if (N[i * mesh.N.size() + n]) {
          NNN.push_back(n);
          break;
        }
      }
    }

    for (int i = 0; i < boxx.size(); ++i) {
      femib::types::dvec<T, d> point = boxx[i];
      femib::types::dvec<T, d> res = {0, 0};
      // if (NNN[i] < mesh.T.size()) {
      femib::types::dtrian<T, d> t = mesh.N[NNN[i]];

      for (int j = 0; j < finite_element.base_functions.size(); ++j) {

        femib::types::F<T, d, d> f = finite_element.base_functions[j];

        std::function<femib::types::dvec<T, d>(femib::types::dvec<T, d>)> g =
            [&](const femib::types::dvec<T, d> &x) {
              return xx(nodes.get_index(j, NNN[i])) *
                     f.x(femib::affine::affine_inv(t, x));
            };
        res += g(point);
      }
      //}

      std::vector<std::vector<float>> uuuu = {{point(0), point(1)},
                                              {res(0), res(1)}};
      uuu.emplace_back(uuuu);
    }
    return uuu;
  }
};
} // namespace femib::finite_element_space
#endif
