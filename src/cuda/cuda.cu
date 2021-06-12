#include <iostream>

#include "../affine/affine.hpp"
#include "../mesh/mesh.hpp"
#include "../types/types.hpp"
#include "cuda.h"
#include "mini-book.h"
#include "spdlog/spdlog.h"

#include <cmath>

#pragma diag_suppress 2527
#pragma diag_suppress 2529
#pragma diag_suppress 2651
#pragma diag_suppress 2653
#pragma diag_suppress 2668
#pragma diag_suppress 2669
#pragma diag_suppress 2670
#pragma diag_suppress 2671
#pragma diag_suppress 2735
#pragma diag_suppress 2737
#pragma diag_suppress 20014

const float EPSILON = std::numeric_limits<float>::epsilon();

void femib::cuda::printSize() {
  spdlog::set_pattern("[%Y-%m-%dT%T] [%l] [%@@%!] %v");
  SPDLOG_INFO("[CUDA stack size] found to be {}", getStackSize());
  SPDLOG_INFO("[CUDA heap  size] found to be {}", getHeapSize());
}

int femib::cuda::getStackSize() {
  size_t size_stack;
  cudaDeviceGetLimit(&size_stack, cudaLimitStackSize);
  return (int)size_stack;
}

int femib::cuda::getHeapSize() {
  size_t size_heap;
  cudaDeviceGetLimit(&size_heap, cudaLimitMallocHeapSize);
  return (int)size_heap;
}

void femib::cuda::setStackSize(int stackSize) {
  cudaDeviceSetLimit(cudaLimitStackSize, stackSize);
}

void femib::cuda::setHeapSize(int heapSize) {
  cudaDeviceSetLimit(cudaLimitMallocHeapSize, heapSize * sizeof(double));
}

template <typename T> T *femib::cuda::copyToDevice(T x) {
  T *X;
  HANDLE_ERROR(cudaMalloc((void **)&X, sizeof(T)));
  HANDLE_ERROR(cudaMemcpy(X, &x, sizeof(T), cudaMemcpyHostToDevice));
  return X;
}

template <typename T> T *femib::cuda::copyVectorToDevice(std::vector<T> x) {
  T *X;
  HANDLE_ERROR(cudaMalloc((void **)&X, sizeof(T) * x.size()));
  HANDLE_ERROR(
      cudaMemcpy(X, x.data(), sizeof(T) * x.size(), cudaMemcpyHostToDevice));
  return X;
}

template <typename T> T femib::cuda::copyToHost(T *X) {
  T x;
  HANDLE_ERROR(cudaMemcpy(&x, X, sizeof(T), cudaMemcpyDeviceToHost));
  return x;
}

#pragma diag_suppress 522
template <typename f, int d>
__host__ __device__ bool
femib::cuda::in_box(const femib::types::dvec<f, d> &P,
                    const femib::types::dtrian<f, d> &T) {
  femib::types::mesh<f, d> mesh = {T};
  femib::types::box<f, d> box = femib::mesh::find_box<f, d>(mesh);
  bool e = true;
  for (int i = 0; e && i < P.size(); ++i) {
    e = e && P(i) > (box[0](i) - EPSILON) && P(i) < (box[1](i) + EPSILON);
  }
  return e;
}

template <typename f, int d>
__host__ __device__ bool
femib::cuda::in_triangle(const femib::types::dvec<f, d> &P,
                         const femib::types::dtrian<f, d> &T) {
  femib::types::dvec<f, d> x = femib::affine::affine_inv<f, d>(T, P);
  return (x(0) >= 0) && (x(1) >= 0) && (x(0) + x(1) <= 1);
}

template <typename f, int d>
__host__ __device__ f distance_point_segment(
    const femib::types::dvec<f, d> &P, const femib::types::dtrian<f, d> &T) {
  femib::types::dvec<f, d> D = T[1] - T[0];
  femib::types::dvec<f, d> E = P - T[0];
  femib::types::dvec<f, d> F = P - T[0];
  f P1P2 = (D.transpose() * D);
  f PP = (E.transpose() * D);
  f dd = PP / P1P2;
  if (dd < 0) {
    return (E.transpose() * E);
  } else if (dd >= 0 && dd <= 1) {
    return ((-E).transpose() * (-E)) - dd * dd * P1P2;
  } else {
    return F.transpose() * F;
  }
}

template <typename f, int d>
__host__ __device__ bool
femib::cuda::accurate(const femib::types::dvec<f, d> &P,
                      const femib::types::dtrian<f, d> &T) {
  if (not femib::cuda::in_box(P, T)) {
    return false;
  }
  if (femib::cuda::in_triangle(P, T)) {
    return true;
  }
  if (false) {
    return false;
  } else if (distance_point_segment(P, {T[0], T[1]}) <= EPSILON) {
    return true;
  } else if (distance_point_segment(P, {T[1], T[2]}) <= EPSILON) {
    return true;
  } else if (distance_point_segment(P, {T[2], T[0]}) <= EPSILON) {
    return true;
  } else {
    return false;
  }
}

template <typename f, int d_>
__host__ __device__ bool
femib::cuda::accurate(const femib::types::dvec<f, d_> &x,
                      femib::types::dvec<f, d_> *t) {
  f a = t[1](0) - t[0](0);
  f b = t[2](0) - t[0](0);
  f c = t[1](1) - t[0](1);
  f d = t[2](1) - t[0](1);

  f X = x(0) - t[0](0);
  f Y = x(1) - t[0](1);

  f det = 1 / (a * d - b * c);

  f x_ = det * (d * X - b * Y);
  f y_ = det * (-c * X + a * Y);

  return (x_ >= 0) && (y_ >= 0) && (x_ + y_ <= 1);
}

template <typename f, int d>
__global__ void parallel_accurate_kernel(femib::types::dvec<f, d> *T,
                                         femib::types::dvec<f, d> *X, bool *N) {
  int blockId = blockIdx.x;
  int threadId = blockId * blockDim.x + threadIdx.x;

  femib::types::dvec<f, d> t[3] = {T[threadId], T[threadId + 1],
                                   T[threadId + 2]};
  N[threadId] = femib::cuda::accurate(X[blockId], t);
}

template <typename f, int d>
__global__ void parallel_accurate_kernel(femib::types::dtrian<f, d> *T,
                                         femib::types::dvec<f, d> *X, bool *N) {
  int blockId = blockIdx.x;
  int threadId = blockId * blockDim.x + threadIdx.x;

  // femib::types::dvec<f, d> t[3] = {T[threadId], T[threadId + 1],
  //                                 T[threadId + 2]};
  N[threadId] = true; // femib::cuda::accurate(X[blockId], t);
}

template <typename f, int d>
__host__ bool
femib::cuda::parallel_accurate(const femib::types::dvec<f, d> &X,
                               const femib::types::dtrian<f, d> &T) {
  femib::types::dvec<f, d> *devT =
      copyVectorToDevice<femib::types::dvec<f, d>>(T);
  femib::types::dvec<f, d> *devX = copyToDevice<femib::types::dvec<f, d>>(X);
  bool N = false;
  bool *devN = copyToDevice(N);
  parallel_accurate_kernel<f, d><<<1, 1>>>(devT, devX, devN);
  bool NN = copyToHost(devN);
  return NN;
}

template <typename f, int d>
__host__ bool serial_accurate(femib::types::dtrian<f, d> T,
                              femib::types::dvec<f, d> X) {

  femib::types::dtrian<f, d> t = T;
  femib::types::dvec<f, d> p = X;
  return femib::cuda::accurate(p, t);
}

template <typename f, int d>
__host__ void femib::cuda::serial_accurate(femib::types::dvec<f, d> *X,
                                           int size_X,
                                           femib::types::dtrian<f, d> *T,
                                           int size_T, bool *N) {
  for (int i = 0; i < size_T; ++i) {
    for (int j = 0; j < size_X; ++j) {
      femib::types::dtrian<f, d> t = T[i];
      femib::types::dvec<f, d> p = X[j];
      N[i * size_T + j] = femib::cuda::accurate(p, t);
    }
  }
}

template <typename f, int d>
__host__ void femib::cuda::parallel_accurate(femib::types::dvec<f, d> *X,
                                             int size_X,
                                             femib::types::dtrian<f, d> *T,
                                             int size_T, bool *N) {

  parallel_accurate_kernel<f, d><<<size_X, size_T>>>(T, X, N);
}

template double *femib::cuda::copyToDevice<double>(double x);
template double femib::cuda::copyToHost<double>(double *x);

template __host__ __device__ bool
femib::cuda::in_box<float, 2>(const femib::types::dvec<float, 2> &P,
                              const femib::types::dtrian<float, 2> &T);
template __host__ __device__ bool
femib::cuda::in_triangle<float, 2>(const femib::types::dvec<float, 2> &P,
                                   const femib::types::dtrian<float, 2> &T);
template __host__ __device__ bool
femib::cuda::accurate<float, 2>(const femib::types::dvec<float, 2> &P,
                                const femib::types::dtrian<float, 2> &T);
template __host__ bool femib::cuda::parallel_accurate<float, 2>(
    const femib::types::dvec<float, 2> &X,
    const femib::types::dtrian<float, 2> &T);

template __host__ void femib::cuda::serial_accurate<float, 2>(
    femib::types::dvec<float, 2> *X, int size_X,
    femib::types::dtrian<float, 2> *T, int size_T, bool *N);
template __host__ void femib::cuda::parallel_accurate<float, 2>(
    femib::types::dvec<float, 2> *X, int size_X,
    femib::types::dtrian<float, 2> *T, int size_T, bool *N);

template <typename T> T *femib::cuda::copyToDevice(T *x, int size) {
  T *X;
  HANDLE_ERROR(cudaMalloc((void **)&X, sizeof(T) * size));
  HANDLE_ERROR(cudaMemcpy(X, x, sizeof(T) * size, cudaMemcpyHostToDevice));
  return X;
}

/*
template <typename T> T *femib::cuda::copyVectorToDevice(std::vector<T> x, int
size) { T *X; HANDLE_ERROR(cudaMalloc((void **)&X, sizeof(T) * x.size()));
  HANDLE_ERROR(
      cudaMemcpy(X, x.data(), sizeof(T) * x.size(), cudaMemcpyHostToDevice));
  return X;
}*/

template <typename T> T *femib::cuda::copyToHost(T *X, int size) {
  T *x;
  HANDLE_ERROR(cudaMemcpy(&x, X, sizeof(T) * size, cudaMemcpyDeviceToHost));
  return x;
}

template femib::types::dvec<float, 2> *
femib::cuda::copyToDevice<femib::types::dvec<float, 2>>(
    femib::types::dvec<float, 2> *x, int size);
template femib::types::dtrian<float, 2> *
femib::cuda::copyToDevice<femib::types::dtrian<float, 2>>(
    femib::types::dtrian<float, 2> *x, int size); // TODO
template bool *femib::cuda::copyToDevice<bool>(bool *x, int size);
template bool *femib::cuda::copyToHost<bool>(bool *x, int size);
