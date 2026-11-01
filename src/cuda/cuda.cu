#include <iostream>

#include "../affine/affine.hpp"
#include "../mesh/mesh.hpp"
#include "../types/types.hpp"
#include "cuda.h"
#include "mini-book.h"
#include "spdlog/spdlog.h"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

void femib::cuda::print_size() {
  SPDLOG_INFO("[CUDA stack size] found to be {}", get_stack_size());
  SPDLOG_INFO("[CUDA heap  size] found to be {}", get_heap_size());
}

size_t femib::cuda::get_stack_size() {
  size_t size_stack;
  cudaDeviceGetLimit(&size_stack, cudaLimitStackSize);
  return size_stack;
}

size_t femib::cuda::get_heap_size() {
  size_t size_heap;
  cudaDeviceGetLimit(&size_heap, cudaLimitMallocHeapSize);
  return size_heap;
}

void femib::cuda::set_stack_size(const size_t stack_size) {
  cudaDeviceSetLimit(cudaLimitStackSize, stack_size);
}

// TODO why mess with this?
void femib::cuda::set_heap_size(const size_t heap_size) {
  cudaDeviceSetLimit(cudaLimitMallocHeapSize, heap_size * sizeof(double));
}

template <typename T> T *femib::cuda::copy_to_device(T *x, const size_t size) {
  T *X;
  HANDLE_ERROR(cudaMalloc((void **)&X, sizeof(T) * size));
  HANDLE_ERROR(cudaMemcpy(X, x, sizeof(T) * size, cudaMemcpyHostToDevice));
  return X;
}

template <typename T> T *femib::cuda::copy_to_host(T *X, const size_t size) {
  T *x = new T[size];
  HANDLE_ERROR(cudaMemcpy(x, X, sizeof(T) * size, cudaMemcpyDeviceToHost));
  return x;
}

template <typename T> void femib::cuda::free_device(T *X) {
  HANDLE_ERROR(cudaFree(X));
}

template <typename T> T *femib::cuda::alloc_device(const size_t size) {
  T *X;
  HANDLE_ERROR(cudaMalloc((void **)&X, sizeof(T) * size));
  return X;
}

template <typename T> void femib::cuda::free_device_quiet(T *X) {
  if (const cudaError_t err = cudaFree(X);
      err != cudaSuccess && err != cudaErrorCudartUnloading) {
    HandleError(err, __FILE__, __LINE__);
  }
}

/******************************************************************************/

template <typename f> __host__ __device__ static f location_tolerance() {
  return std::sqrt(std::numeric_limits<f>::epsilon()) / 4;
}

template <typename f, int d>
__host__ __device__ static bool in_box_(const femib::types::dvec<f, d> &P,
                                        const femib::types::dvec<f, d> *T) {
  for (int i = 0; i < d; ++i) {
    f min = T[0](i), max = T[0](i);
    for (int k = 1; k < d + 1; ++k) {
      if (T[k](i) < min)
        min = T[k](i);
      if (T[k](i) > max)
        max = T[k](i);
    }
    // wide enough for every point accurate_ accepts
    const f slack = 2 * location_tolerance<f>() * (max - min);
    if (!(P(i) >= min - slack && P(i) <= max + slack)) {
      return false;
    }
  }
  return true;
}

template <typename f, int d>
__host__ __device__ static bool accurate_(const femib::types::dvec<f, d> &x,
                                          const femib::types::dvec<f, d> *t) {
  if (!in_box_<f, d>(x, t)) {
    return false;
  }
  f a = t[1](0) - t[0](0);
  f b = t[2](0) - t[0](0);
  f c = t[1](1) - t[0](1);
  f dd_ = t[2](1) - t[0](1);
  f X = x(0) - t[0](0);
  f Y = x(1) - t[0](1);
  f det = 1 / (a * dd_ - b * c);
  f x_ = det * (dd_ * X - b * Y);
  f y_ = det * (-c * X + a * Y);
  const f tol = location_tolerance<f>();
  return x_ >= -tol && y_ >= -tol && x_ + y_ <= 1 + tol;
}

template <typename f, int d>
__host__ bool femib::cuda::in_box(const femib::types::dvec<f, d> &P,
                                  const femib::types::dtrian<f, d> &T) {
  return in_box_<f, d>(P, T.data());
}

template <typename f, int d>
__host__ bool femib::cuda::in_triangle(const femib::types::dvec<f, d> &P,
                                       const femib::types::dtrian<f, d> &T) {
  femib::types::dvec<f, d> x = femib::affine::affine_inv<f, d>(T, P);
  return (x(0) >= 0) && (x(1) >= 0) && (x(0) + x(1) <= 1);
}

template <typename f, int d>
__host__ bool femib::cuda::accurate(const femib::types::dvec<f, d> &P,
                                    const femib::types::dtrian<f, d> &T) {
  return accurate_<f, d>(P, T.data());
}

template <typename f, int d>
__host__ void femib::cuda::serial_accurate(femib::types::dvec<f, d> *X,
                                           const size_t size_X,
                                           femib::types::dtrian<f, d> *T,
                                           const size_t size_T, bool *N) {
  for (size_t j = 0; j < size_X; ++j) {
    for (size_t i = 0; i < size_T; ++i) {
      femib::types::dtrian<f, d> t = T[i];
      femib::types::dvec<f, d> p = X[j];
      N[j * size_T + i] = femib::cuda::accurate(p, t);
    }
  }
}

template <typename f, int d>
__global__ static void
parallel_accurate_kernel(femib::types::dtrian_<f, d> *T, const size_t size_T,
                         femib::types::dvec<f, d> *X, bool *N) {
  const size_t point_idx = blockIdx.x;
  const size_t tri_idx = size_t{blockIdx.y} * blockDim.x + threadIdx.x;
  if (tri_idx >= size_T) {
    return;
  }
  femib::types::dvec<f, d> p = X[point_idx];
  const bool n = accurate_<f, d>(p, T[tri_idx]);
  N[point_idx * size_T + tri_idx] = n;
}

template <typename f, int d>
__host__ void femib::cuda::parallel_accurate(femib::types::dvec<f, d> *X,
                                             const size_t size_X,
                                             femib::types::dtrian_<f, d> *T,
                                             const size_t size_T, bool *N) {
  constexpr unsigned int threads_per_block = 256; //  hardware limit = 1024
  const size_t blocks_per_point =
      (size_T + threads_per_block - 1) / threads_per_block;

  int device;
  HANDLE_ERROR(cudaGetDevice(&device));
  int max_grid_x;
  int max_grid_y;
  HANDLE_ERROR(
      cudaDeviceGetAttribute(&max_grid_x, cudaDevAttrMaxGridDimX, device));
  HANDLE_ERROR(
      cudaDeviceGetAttribute(&max_grid_y, cudaDevAttrMaxGridDimY, device));
  if (size_X > static_cast<size_t>(max_grid_x)) {
    throw std::invalid_argument(
        "femib::cuda::parallel_accurate: " + std::to_string(size_X) +
        " points exceed the device's grid limit of " +
        std::to_string(max_grid_x));
  }
  if (blocks_per_point > static_cast<size_t>(max_grid_y)) {
    throw std::invalid_argument(
        "femib::cuda::parallel_accurate: " + std::to_string(size_T) +
        " triangles exceed the device's grid limit of " +
        std::to_string(static_cast<size_t>(max_grid_y) * threads_per_block));
  }
  if (size_X == 0 || size_T == 0) {
    return;
  }

  // both fit in unsigned int after the checks above
  const dim3 grid(static_cast<unsigned int>(size_X),
                  static_cast<unsigned int>(blocks_per_point));
  parallel_accurate_kernel<f, d><<<grid, threads_per_block>>>(T, size_T, X, N);
}

/******************************************************************************/

template double *femib::cuda::copy_to_device<double>(double *x, size_t size);
template double *femib::cuda::copy_to_host<double>(double *x, size_t size);
template void femib::cuda::free_device<double>(double *x);

template femib::types::dvec<float, 2> *
femib::cuda::copy_to_device<femib::types::dvec<float, 2>>(
    femib::types::dvec<float, 2> *x, size_t size);
template femib::types::dtrian_<float, 2> *
femib::cuda::copy_to_device<femib::types::dtrian_<float, 2>>(
    femib::types::dtrian_<float, 2> *x, size_t size);
template bool *femib::cuda::copy_to_device<bool>(bool *x, size_t size);
template bool *femib::cuda::copy_to_host<bool>(bool *x, size_t size);
template bool *femib::cuda::alloc_device<bool>(size_t size);

template void femib::cuda::free_device<femib::types::dvec<float, 2>>(
    femib::types::dvec<float, 2> *x);
template void femib::cuda::free_device<femib::types::dtrian_<float, 2>>(
    femib::types::dtrian_<float, 2> *x);
template void femib::cuda::free_device<bool>(bool *x);
template void femib::cuda::free_device_quiet<femib::types::dtrian_<float, 2>>(
    femib::types::dtrian_<float, 2> *x);

/******************************************************************************/

template __host__ bool
femib::cuda::in_box<float, 2>(const femib::types::dvec<float, 2> &P,
                              const femib::types::dtrian<float, 2> &T);
template __host__ bool
femib::cuda::in_triangle<float, 2>(const femib::types::dvec<float, 2> &P,
                                   const femib::types::dtrian<float, 2> &T);
template __host__ bool
femib::cuda::accurate<float, 2>(const femib::types::dvec<float, 2> &P,
                                const femib::types::dtrian<float, 2> &T);
template __host__ void femib::cuda::serial_accurate<float, 2>(
    femib::types::dvec<float, 2> *X, size_t size_X,
    femib::types::dtrian<float, 2> *T, size_t size_T, bool *N);
template __host__ void femib::cuda::parallel_accurate<float, 2>(
    femib::types::dvec<float, 2> *X, size_t size_X,
    femib::types::dtrian_<float, 2> *T, size_t size_T, bool *N);

/******************************************************************************/

template femib::types::dvec<double, 2> *
femib::cuda::copy_to_device<femib::types::dvec<double, 2>>(
    femib::types::dvec<double, 2> *x, size_t size);
template femib::types::dtrian_<double, 2> *
femib::cuda::copy_to_device<femib::types::dtrian_<double, 2>>(
    femib::types::dtrian_<double, 2> *x, size_t size);
template void femib::cuda::free_device<femib::types::dvec<double, 2>>(
    femib::types::dvec<double, 2> *x);
template void femib::cuda::free_device<femib::types::dtrian_<double, 2>>(
    femib::types::dtrian_<double, 2> *x);
template void femib::cuda::free_device_quiet<femib::types::dtrian_<double, 2>>(
    femib::types::dtrian_<double, 2> *x);

template __host__ bool
femib::cuda::in_box<double, 2>(const femib::types::dvec<double, 2> &P,
                               const femib::types::dtrian<double, 2> &T);
template __host__ bool
femib::cuda::in_triangle<double, 2>(const femib::types::dvec<double, 2> &P,
                                    const femib::types::dtrian<double, 2> &T);
template __host__ bool
femib::cuda::accurate<double, 2>(const femib::types::dvec<double, 2> &P,
                                 const femib::types::dtrian<double, 2> &T);
template __host__ void femib::cuda::serial_accurate<double, 2>(
    femib::types::dvec<double, 2> *X, size_t size_X,
    femib::types::dtrian<double, 2> *T, size_t size_T, bool *N);
template __host__ void femib::cuda::parallel_accurate<double, 2>(
    femib::types::dvec<double, 2> *X, size_t size_X,
    femib::types::dtrian_<double, 2> *T, size_t size_T, bool *N);

/******************************************************************************/
