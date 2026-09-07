#include <iostream>

#include "../affine/affine.hpp"
#include "../mesh/mesh.hpp"
#include "../types/types.hpp"
#include "cuda.h"
#include "mini-book.h"
#include "spdlog/spdlog.h"

#include <cmath>
#include <limits>

void femib::cuda::printSize() {
  SPDLOG_INFO("[CUDA stack size] found to be {}", getStackSize());
  SPDLOG_INFO("[CUDA heap  size] found to be {}", getHeapSize());
}

size_t femib::cuda::getStackSize() {
  size_t size_stack;
  cudaDeviceGetLimit(&size_stack, cudaLimitStackSize);
  return size_stack;
}

size_t femib::cuda::getHeapSize() {
  size_t size_heap;
  cudaDeviceGetLimit(&size_heap, cudaLimitMallocHeapSize);
  return size_heap;
}

void femib::cuda::setStackSize(size_t stackSize) {
  cudaDeviceSetLimit(cudaLimitStackSize, stackSize);
}

// TODO why mess with this?
void femib::cuda::setHeapSize(size_t heapSize) {
  cudaDeviceSetLimit(cudaLimitMallocHeapSize, heapSize * sizeof(double));
}

template <typename T> T *femib::cuda::copyToDevice(T *x, int size) {
  T *X;
  HANDLE_ERROR(cudaMalloc((void **)&X, sizeof(T) * size));
  HANDLE_ERROR(cudaMemcpy(X, x, sizeof(T) * size, cudaMemcpyHostToDevice));
  return X;
}

template <typename T> T *femib::cuda::copyToHost(T *X, int size) {
  T *x = new T[size];
  HANDLE_ERROR(cudaMemcpy(x, X, sizeof(T) * size, cudaMemcpyDeviceToHost));
  return x;
}

template <typename T> void femib::cuda::freeDevice(T *X) {
  HANDLE_ERROR(cudaFree(X));
}

/******************************************************************************/

template <typename f, int d>
__host__ bool femib::cuda::in_box(const femib::types::dvec<f, d> &P,
                                  const femib::types::dtrian<f, d> &T) {
  f EPSILON = std::numeric_limits<f>::epsilon();
  femib::types::mesh<f, d> mesh = {T};
  femib::types::box<f, d> box = femib::mesh::find_box<f, d>(mesh);
  bool e = true;
  for (int i = 0; e && i < P.size(); ++i) {
    e = e && P(i) > (box[0](i) - EPSILON) && P(i) < (box[1](i) + EPSILON);
  }
  return e;
}

template <typename f, int d>
__device__ bool in_box_(const femib::types::dvec<f, d> &P,
                        femib::types::dvec<f, d> *T) {
  f EPSILON = std::numeric_limits<f>::epsilon();
  for (int i = 0; i < d; ++i) {
    f min = T[0](i), max = T[0](i);
    for (int k = 1; k < d + 1; ++k) { // TODO d+1?
      if (T[k](i) < min)
        min = T[k](i);
      if (T[k](i) > max)
        max = T[k](i);
    }
    if (!(P(i) > min - EPSILON && P(i) < max + EPSILON)) {
      return false;
    }
  }
  return true;
}

template <typename f, int d>
__host__ bool femib::cuda::in_triangle(const femib::types::dvec<f, d> &P,
                                       const femib::types::dtrian<f, d> &T) {
  femib::types::dvec<f, d> x = femib::affine::affine_inv<f, d>(T, P);
  return (x(0) >= 0) && (x(1) >= 0) && (x(0) + x(1) <= 1);
}

template <typename f, int d>
__host__ f distance_point_segment(const femib::types::dvec<f, d> &P,
                                  const femib::types::dtrian<f, d> &T) {
  femib::types::dvec<f, d> D = T[1] - T[0];
  femib::types::dvec<f, d> E = P - T[0];
  femib::types::dvec<f, d> F = P - T[1];
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
__device__ f distance_point_segment_(const femib::types::dvec<f, d> &P,
                                     const femib::types::dvec<f, d> &A,
                                     const femib::types::dvec<f, d> &B) {
  femib::types::dvec<f, d> D = B - A;
  femib::types::dvec<f, d> E = P - A;
  f P1P2 = D.dot(D);
  f dd = E.dot(D) / P1P2;
  if (dd < 0) {
    return E.dot(E);
  } else if (dd <= 1) {
    return E.dot(E) - dd * dd * P1P2;
  } else {
    femib::types::dvec<f, d> F = P - B;
    return F.dot(F);
  }
}

// TODO what algoritm is this? source?
template <typename f, int d>
__host__ bool femib::cuda::accurate(const femib::types::dvec<f, d> &P,
                                    const femib::types::dtrian<f, d> &T) {
  f EPSILON = std::numeric_limits<f>::epsilon();
  if (not femib::cuda::in_box(P, T)) {
    return false;
  }
  if (femib::cuda::in_triangle(P, T)) {
    return true;
  }
  // TODO eh
  if (false) {
    return false;
  } else if (distance_point_segment(P, {T[0], T[1]}) <= EPSILON * EPSILON) {
    return true;
  } else if (distance_point_segment(P, {T[1], T[2]}) <= EPSILON * EPSILON) {
    return true;
  } else if (distance_point_segment(P, {T[2], T[0]}) <= EPSILON * EPSILON) {
    return true;
  } else {
    return false;
  }
}

template <typename f, int d>
__device__ bool accurate_(const femib::types::dvec<f, d> &x,
                          femib::types::dvec<f, d> *t) {
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
  if (x_ >= 0 && y_ >= 0 && x_ + y_ <= 1) {
    return true;
  }
  f EPSILON = std::numeric_limits<f>::epsilon();
  if (distance_point_segment_<f, d>(x, t[0], t[1]) <= EPSILON * EPSILON)
    return true;
  if (distance_point_segment_<f, d>(x, t[1], t[2]) <= EPSILON * EPSILON)
    return true;
  if (distance_point_segment_<f, d>(x, t[2], t[0]) <= EPSILON * EPSILON)
    return true;
  return false;
}

/******************************************************************************/

template <typename f, int d>
__host__ void femib::cuda::serial_accurate(femib::types::dvec<f, d> *X,
                                           int size_X,
                                           femib::types::dtrian<f, d> *T,
                                           int size_T, bool *N) {
  for (int j = 0; j < size_X; ++j) {
    for (int i = 0; i < size_T; ++i) {
      femib::types::dtrian<f, d> t = T[i];
      femib::types::dvec<f, d> p = X[j];
      N[j * size_T + i] = femib::cuda::accurate(p, t);
    }
  }
}

template <typename f, int d>
__global__ void parallel_accurate_kernel(femib::types::dtrian_<f, d> *T,
                                         femib::types::dvec<f, d> *X, bool *N) {
  int blockId = blockIdx.x;
  int threadId = blockId * blockDim.x + threadIdx.x;
  femib::types::dvec<f, d> p = X[blockId];
  bool n = accurate_<f, d>(p, T[threadIdx.x]);
  N[threadId] = n;
}

template <typename f, int d>
__host__ void femib::cuda::parallel_accurate(femib::types::dvec<f, d> *X,
                                             int size_X,
                                             femib::types::dtrian_<f, d> *T,
                                             int size_T, bool *N) {
  parallel_accurate_kernel<f, d><<<size_X, size_T>>>(T, X, N);
}

/******************************************************************************/

template double *femib::cuda::copyToDevice<double>(double *x, int size);
template double *femib::cuda::copyToHost<double>(double *x, int size);

template femib::types::dvec<float, 2> *
femib::cuda::copyToDevice<femib::types::dvec<float, 2>>(
    femib::types::dvec<float, 2> *x, int size);
template femib::types::dtrian_<float, 2> *
femib::cuda::copyToDevice<femib::types::dtrian_<float, 2>>(
    femib::types::dtrian_<float, 2> *x, int size);
template bool *femib::cuda::copyToDevice<bool>(bool *x, int size);
template bool *femib::cuda::copyToHost<bool>(bool *x, int size);

template void femib::cuda::freeDevice<femib::types::dvec<float, 2>>(
    femib::types::dvec<float, 2> *x);
template void femib::cuda::freeDevice<femib::types::dtrian_<float, 2>>(
    femib::types::dtrian_<float, 2> *x);
template void femib::cuda::freeDevice<bool>(bool *x);

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
    femib::types::dvec<float, 2> *X, int size_X,
    femib::types::dtrian<float, 2> *T, int size_T, bool *N);
template __host__ void femib::cuda::parallel_accurate<float, 2>(
    femib::types::dvec<float, 2> *X, int size_X,
    femib::types::dtrian_<float, 2> *T, int size_T, bool *N);

/******************************************************************************/
