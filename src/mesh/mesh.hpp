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
  size_t n_triangles = mesh.N.size();

  femib::types::dtrian_<T, d> *device_triangles;
  if (mesh.device_triangles_cache) {
    device_triangles = static_cast<femib::types::dtrian_<T, d> *>(
        mesh.device_triangles_cache.get());
  } else {
    femib::types::dtrian_<T, d> *host_triangles =
        femib::types::vector_dtrian2pointer_dtrian_<T, d>(mesh.N);
    device_triangles = femib::cuda::copy_to_device<femib::types::dtrian_<T, d>>(
        host_triangles, n_triangles);
    delete[] host_triangles;
    mesh.device_triangles_cache =
        std::shared_ptr<void>(device_triangles, [](void *p) {
          femib::cuda::free_device_quiet<femib::types::dtrian_<T, d>>(
              static_cast<femib::types::dtrian_<T, d> *>(p));
        });
  }
  femib::types::dvec<T, d> *device_points =
      femib::cuda::copy_to_device<femib::types::dvec<T, d>>(points.data(),
                                                            points.size());

  bool *device_hits =
      femib::cuda::alloc_device<bool>(points.size() * n_triangles);

  femib::cuda::parallel_accurate<T, d>(
      device_points, points.size(), device_triangles, n_triangles, device_hits);
  const bool *hits =
      femib::cuda::copy_to_host<bool>(device_hits, points.size() * n_triangles);

  std::vector<int> containing_triangles;
  containing_triangles.reserve(points.size());

  for (size_t i = 0; i < points.size(); ++i) {
    int found = -1;
    for (int n = 0; n < std::ssize(mesh.N); ++n) {
      if (hits[i * n_triangles + n]) {
        found = n;
        break;
      }
    }
    containing_triangles.push_back(found);
  }
  cuda::free_device<femib::types::dvec<T, d>>(device_points);
  cuda::free_device<bool>(device_hits);
  delete[] hits;
  hits = nullptr;
  return containing_triangles;
}

template <typename T, int d>
T integrate(const femib::gauss::rule<T, d> &rule,
            const std::function<T(femib::types::dvec<T, d>)> &f,
            const femib::types::dtrian<T, d> &t) {
  std::function<T(femib::types::dvec<T, d>)> g =
      [&t, &f](const femib::types::dvec<T, d> &x) {
        return femib::affine::abs_jacobian_determinant(t) *
               f(femib::affine::reference_to_physical(t, x));
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
  return std::transform_reduce(std::execution::par_unseq, mesh.N.begin(),
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
femib::types::box<T, d>
find_box(const std::vector<femib::types::dvec<T, d>> &points) {
  if (points.empty()) {
    throw std::invalid_argument("find_box: mesh is empty");
  }
  femib::types::dvec<T, d> b1 = points[0];
  femib::types::dvec<T, d> b2 = points[0];
  for (size_t n = 1; n < points.size(); ++n) {
    b1 = b1.array().min(points[n].array());
    b2 = b2.array().max(points[n].array());
  }
  return femib::types::box<T, d>{.bottom = b1, .top = b2};
}

template <typename T, int d>
femib::types::box<T, d> find_box(const femib::types::mesh<T, d> &m) {
  return find_box<T, d>(m.P);
}

template <typename T, int d>
std::vector<femib::types::dvec<T, d>>
build_uniform_grid(const femib::types::box<T, d> &box, T delta) {
  T x_min = box.bottom(0);
  T x_max = box.top(0);
  T y_min = box.bottom(1);
  T y_max = box.top(1);
  std::vector<femib::types::dvec<T, d>> grid;
  for (T x = x_min; x <= x_max; x += delta) {
    for (T y = y_min; y <= y_max; y += delta) {
      femib::types::dvec<T, d> w = {x, y};
      grid.emplace_back(w);
    }
  }
  return grid;
}

} // namespace femib::mesh
#endif
