#ifndef MESH_HPP_INCLUDED_
#define MESH_HPP_INCLUDED_

#include "../affine/affine.hpp"
#include "../cuda/cuda.h"
#include "../gauss/gauss.hpp"
#include "../read/read.hpp"
#include "../types/types.hpp"
#include <execution>
#include <functional>
#include <numeric>
#include <stdexcept>
#include <string>

namespace femib::mesh {

template <typename T, int d>
std::vector<int> find_points(const femib::types::mesh<T, d> &mesh,
                             std::vector<types::dvec<T, d>> &points) {
  size_t meshSize = mesh.N.size();

  // TODO we can load this during init, only points changes, the mesh is always
  //  the same
  femib::types::dtrian_<T, d> *T_ =
      femib::types::vector_dtrian2pointer_dtrian_<T, d>(mesh.N);
  femib::types::dtrian_<T, d> *devT =
      femib::cuda::copyToDevice<femib::types::dtrian_<T, d>>(T_, meshSize);
  femib::types::dvec<T, d> *devX =
      femib::cuda::copyToDevice<femib::types::dvec<T, d>>(points.data(),
                                                          points.size());
  std::unique_ptr<bool[]> Nn(
      new bool[points.size() *
               meshSize]); // was: bool Nn[points.size() * meshSize];

  bool *devN =
      femib::cuda::copyToDevice<bool>(Nn.get(), points.size() * meshSize);

  femib::cuda::parallel_accurate<T, d>(devX, points.size(), devT, meshSize,
                                       devN);
  bool *NN = femib::cuda::copyToHost<bool>(devN, points.size() * meshSize);

  std::vector<int> NNN;
  NNN.reserve(points.size() * meshSize);

  for (int i = 0; i < points.size(); ++i) {
    int found = -1;
    for (int n = 0; n < meshSize; ++n) {
      if (NN[i * meshSize + n]) {
        found = n;
        break;
      }
    }
    NNN.push_back(found);
  }
  cuda::freeDevice<femib::types::dtrian_<T, d>>(devT);
  cuda::freeDevice<femib::types::dvec<T, d>>(devX);
  cuda::freeDevice<bool>(devN);
  delete[] T_;
  T_ = nullptr;
  delete[] NN;
  NN = nullptr;
  return NNN;
}

template <typename T, int d>
T integrate(const femib::gauss::rule<T, d> &rule,
            const std::function<T(femib::types::dvec<T, d>)> &f,
            const femib::types::dtrian<T, d> &t) {
  std::function<T(femib::types::dvec<T, d>)> g =
      [&t, &f](const femib::types::dvec<T, d> &x) {
        return femib::affine::affineBdet(t) * f(femib::affine::affine(t, x));
      };
  return femib::gauss::integrate<T, d>(rule, g);
}

template <typename T, int d>
T integrate(const femib::gauss::rule<T, d> &rule,
            const std::function<T(femib::types::dvec<T, d>)> &f,
            const femib::types::mesh<T, d> &mesh) {
  auto unary_op = [&rule, &f](const femib::types::dtrian<T, d> &t) {
    return femib::mesh::integrate(rule, f, t);
  };
  return std::transform_reduce(std::execution::seq, mesh.N.begin(),
                               mesh.N.end(), T(0.0), std::plus<>(), unary_op);
}

template <typename T, int d>
femib::types::mesh<T, d> read(const std::string &filename_p,
                              const std::string &filename_t,
                              const std::string &filename_e) {
  femib::types::mesh<T, d> mesh =
      read_mesh_file<T, d>(filename_p, filename_t, filename_e);
  return mesh;
}

template <typename T, int d>
femib::types::box<T, d> find_box(const femib::types::mesh<T, d> &m) {
  if (m.P.empty()) {
    throw std::invalid_argument("find_box: mesh is empty");
  }
  femib::types::box<T, d> box;
  femib::types::dvec<T, d> b1 = m.P[0];
  femib::types::dvec<T, d> b2 = m.P[0];
  for (int n = 1; n < m.P.size(); ++n) {
    b1 = b1.array().min(m.P[n].array());
    b2 = b2.array().max(m.P[n].array());
  }
  box.emplace_back(b1);
  box.emplace_back(b2);
  return box;
}

template <typename T, int d>
femib::types::box<T, d> lin_spaced(
    const femib::types::box<T, d> &b,
    T delta) { // TODO build_uniform_grid, also a box is exactly 2 points (top,
               // bottom), a grid is just a bunch of points, so std::vector
  T x_min = b[0](0);
  T x_max = b[1](0);
  T y_min = b[0](1);
  T y_max = b[1](1);
  femib::types::box<T, d> box;
  for (T x = x_min; x <= x_max; x += delta) {
    for (T y = y_min; y <= y_max; y += delta) {
      femib::types::dvec<T, d> w = {x, y};
      box.emplace_back(w);
    }
  }
  return box;
}

} // namespace femib::mesh
#endif
