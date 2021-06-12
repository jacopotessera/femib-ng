#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../src/types/types.hpp"
#include "cuda.h"
#include <doctest/doctest.h>

femib::types::dtrian<float, 2> T = {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}};
femib::types::dvec<float, 2> P1 = {0.5, 0.5};
femib::types::dvec<float, 2> P2 = {1.5, 0.5};
femib::types::dvec<float, 2> P3 = {0.0, 0.0};
femib::types::dvec<float, 2> P4 = {-0.5, -0.5};
femib::types::dvec<float, 2> P5 = {0.7, 0.7};

femib::types::dtrian<float, 2> Ts[] = {T};
femib::types::dvec<float, 2> Ps[] = {P1, P2, P3, P4, P5};

TEST_CASE("testing cuda size") {
  femib::cuda::setStackSize(FEMIB_CUDA_STACK_SIZE);
  femib::cuda::setHeapSize(FEMIB_CUDA_HEAP_SIZE);
  femib::cuda::printSize();
  WARN(femib::cuda::getStackSize() == FEMIB_CUDA_STACK_SIZE);
  WARN(femib::cuda::getHeapSize() == FEMIB_CUDA_HEAP_SIZE);
}

TEST_CASE("testing cuda copy") {
  double x = 10;
  double *X = femib::cuda::copyToDevice<double>(x);
  double y = femib::cuda::copyToHost<double>(X);
  CHECK(x == y);
}

TEST_CASE("testing cuda in_box") {
  CHECK(femib::cuda::in_box<float, 2>(P1, T));
  CHECK_FALSE(femib::cuda::in_box<float, 2>(P2, T));
  CHECK(femib::cuda::in_box<float, 2>(P3, T));
  CHECK_FALSE(femib::cuda::in_box<float, 2>(P4, T));
  CHECK(femib::cuda::in_box<float, 2>(P5, T));
}

TEST_CASE("testing cuda in_triangle") {
  CHECK(femib::cuda::in_triangle<float, 2>(P1, T));
  CHECK_FALSE(femib::cuda::in_triangle<float, 2>(P2, T));
  CHECK(femib::cuda::in_triangle<float, 2>(P3, T));
  CHECK_FALSE(femib::cuda::in_triangle<float, 2>(P4, T));
  CHECK_FALSE(femib::cuda::in_triangle<float, 2>(P5, T));
}

TEST_CASE("testing cuda accurate") {
  CHECK(femib::cuda::accurate<float, 2>(P1, T));
  CHECK_FALSE(femib::cuda::accurate<float, 2>(P2, T));
  CHECK(femib::cuda::accurate<float, 2>(P3, T));
  CHECK_FALSE(femib::cuda::accurate<float, 2>(P4, T));
  CHECK_FALSE(femib::cuda::accurate<float, 2>(P5, T));
}

TEST_CASE("testing cuda serial_accurate") {
  bool N[5];
  femib::cuda::serial_accurate<float, 2>(Ps, 5, Ts, 1, N);
  CHECK(N[0]);
  CHECK_FALSE(N[1]);
  CHECK(N[2]);
  CHECK_FALSE(N[3]);
  CHECK_FALSE(N[4]);
}

TEST_CASE("testing cuda parallel_accurate") {
  bool N[5];
  femib::types::dtrian<float, 2> *devT =
      femib::cuda::copyToDevice<femib::types::dtrian<float, 2>>(Ts, 1);
  femib::types::dvec<float, 2> *devX =
      femib::cuda::copyToDevice<femib::types::dvec<float, 2>>(Ps, 5);
  bool *devN = femib::cuda::copyToDevice<bool>(N, 5);
  femib::cuda::parallel_accurate<float, 2>(devX, 5, devT, 1, devN);
  bool *NN = femib::cuda::copyToHost<bool>(devN, 5);
  CHECK(N[0]);
  CHECK_FALSE(N[1]);
  CHECK(N[2]);
  CHECK_FALSE(N[3]);
  CHECK_FALSE(N[4]);
}
