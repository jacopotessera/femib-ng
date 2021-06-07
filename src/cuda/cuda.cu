#include <iostream>

#include "../mesh/mesh.hpp"
#include "../types/types.hpp"
#include "../affine/affine.hpp"
#include "cuda.h"
#include "mini-book.h"
#include "spdlog/spdlog.h"

#include <cmath>

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

template <typename T> T femib::cuda::copyToHost(T *X) {
  T x;
  HANDLE_ERROR(cudaMemcpy(&x, X, sizeof(T), cudaMemcpyDeviceToHost));
  return x;
}

template <typename f, int d>
__host__ __device__ bool
femib::cuda::in_box(const femib::types::dvec<f, d> &P,
                    const femib::types::dtrian<f, d> &T) {
  femib::types::mesh<f, d> mesh = {T};
  femib::types::box<f, d> box = femib::mesh::find_box<f, d>(mesh);
  bool e = true;
  for (int i = 0; e && i < P.size(); ++i) {
    e = e && P(i) > (box[0](i) - sqrtf(EPSILON)) &&
        P(i) < (box[1](i) + sqrtf(EPSILON));
  }
  return e;
}

template <typename f, int d>
__host__ __device__ bool
femib::cuda::in_triangle(const femib::types::dvec<f, d> &P,
                         const femib::types::dtrian<f, d> &T) {

  femib::types::dvec<f, d> b[2];
  b[0] = T[1] - T[0];
  b[1] = T[2] - T[0];
  femib::types::dvec<f, d> p = P - T[0];
  femib::types::dmat<f, d> M;
  for (int i = 0; i < M.rows(); ++i) {
    for (int j = 0; j < M.cols(); ++j) {
      M(i, j) = b[j](i);
    }
  }

  femib::types::dvec<f, d> x = M.inverse() * p;
  return (x(0) >= 0) && (x(1) >= 0) && ((x(0) + x(1)) <= 1);
}

/*
__host__ __device__ template <typename f, int d>
bool accurate(const femib::types::dvec<f, d> &P,
              const femib::types::dtrian<f, d> &T) {
  bool N;
  if (not in_box(P, T)) {
    N = 0;
  } else {
    if (in_triangle(P, T)) {
      N = 1;
    } else {
      if (P.size == 1 && T.size == 2) {
        if (false) {
          N = 0;
        } else if (distancePointPoint(P, T(0)) <= EPSILON) {
          N = 1;
        } else if (distancePointPoint(P, T(1)) <= EPSILON) {
          N = 1;
        } else {
          N = 0;
        }
      } else if (P.size == 2 && T.size == 3) {
        if (false) {
          N = 0;
        } else if (distancePointSegment(P, {T(0), T(1)}) <= EPSILON) {
          N = 1;
        } else if (distancePointSegment(P, {T(1), T(2)}) <= EPSILON) {
          N = 1;
        } else if (distancePointSegment(P, {T(2), T(0)}) <= EPSILON) {
          N = 1;
        } else {
          N = 0;
        }
      } else if (P.size == 3 && T.size == 4) {
        if (false) {
          N = 0;
        } else if (distancePointTriangle(P, {T(0), T(1), T(2)}) <= EPSILON) {
          N = 1;
        } else if (distancePointTriangle(P, {T(0), T(1), T(3)}) <= EPSILON) {
          N = 1;
        } else if (distancePointTriangle(P, {T(1), T(2), T(3)}) <= EPSILON) {
          N = 1;
        } else if (distancePointTriangle(P, {T(2), T(0), T(3)}) <= EPSILON) {
          N = 1;
        } else {
          N = 0;
        }
      } else {
        N = 0;
      }
    }
  }
  return N;
}*/

template double *femib::cuda::copyToDevice<double>(double x);
template double femib::cuda::copyToHost<double>(double *x);

template __host__ __device__ bool
femib::cuda::in_box<float, 2>(const femib::types::dvec<float, 2> &P,
                              const femib::types::dtrian<float, 2> &T);
template __host__ __device__ bool
femib::cuda::in_triangle<float, 2>(const femib::types::dvec<float, 2> &P,
                                   const femib::types::dtrian<float, 2> &T);
