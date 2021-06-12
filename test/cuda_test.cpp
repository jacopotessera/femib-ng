#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../src/types/types.hpp"
#include "cuda.h"
#include <doctest/doctest.h>
#include <vector>

femib::types::dtrian<float, 2> T = {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}};
femib::types::dtrian<float, 2> T2 = {{1.0, 0.0}, {1.0, 0.0}, {0.0, -1.0}};
femib::types::dtrian<float, 2> T3 = {{2.0, 0.0}, {-1.0, 0.0}, {0.0, 1.0}};
femib::types::dtrian<float, 2> T4 = {{3.0, 0.0}, {-2.0, 2.0}, {0.0, -1.0}};
femib::types::dtrian<float, 2> T5 = {
    {-4.0, 0.0}, {-1.0, -20.0}, {-10.0, -11.0}};
femib::types::dvec<float, 2> P1 = {0.5, 0.5};
femib::types::dvec<float, 2> P2 = {1.5, 0.5};
femib::types::dvec<float, 2> P3 = {0.0, 0.0};
femib::types::dvec<float, 2> P4 = {-5.0, -10.0};
femib::types::dvec<float, 2> P5 = {0.7, 0.7};

femib::types::dtrian<float, 2> Ts[] = {T};
femib::types::dtrian<float, 2> Tss[] = {T, T2, T3, T4, T5};
femib::types::dvec<float, 2> Ps[] = {P1, P2, P3, P4, P5};
femib::types::dvec<float, 2> Pss[] = {P1};

TEST_CASE("testing cuda size") {
  femib::cuda::setStackSize(FEMIB_CUDA_STACK_SIZE);
  femib::cuda::setHeapSize(FEMIB_CUDA_HEAP_SIZE);
  femib::cuda::printSize();
  WARN(femib::cuda::getStackSize() == FEMIB_CUDA_STACK_SIZE);
  WARN(femib::cuda::getHeapSize() == FEMIB_CUDA_HEAP_SIZE);
}

TEST_CASE("testing cuda copy") {
  double x = 10;
  double *X = femib::cuda::copyToDevice<double>(&x, 1);
  double *y = femib::cuda::copyToHost<double>(X, 1);
  CHECK(x == *y);
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
  bool N[25];
  femib::cuda::serial_accurate<float, 2>(Ps, 5, Tss, 5, N);
  CHECK(N[0]);
  CHECK_FALSE(N[1]);
  CHECK(N[2]);
  CHECK_FALSE(N[3]);
  CHECK_FALSE(N[4]);
  // ...
  CHECK(N[23]);
  CHECK_FALSE(N[24]);
}

TEST_CASE("testing cuda parallel_accurate") {
  bool N[25];

  std::vector<femib::types::dtrian<float, 2>> Tss_ = {T, T2, T3, T4, T5};
  femib::types::dtrian_<float, 2> *Tss__ =
      femib::types::vector_dtrian2pointer_dtrian_<float, 2>(Tss_);

  femib::types::dtrian_<float, 2> *devT =
      femib::cuda::copyToDevice<femib::types::dtrian_<float, 2>>(Tss__, 5);
  femib::types::dvec<float, 2> *devX =
      femib::cuda::copyToDevice<femib::types::dvec<float, 2>>(Ps, 5);
  bool *devN = femib::cuda::copyToDevice<bool>(N, 25);
  femib::cuda::parallel_accurate<float, 2>(devX, 5, devT, 5, devN);
  bool *NN;
  NN = femib::cuda::copyToHost<bool>(devN, 25);
  CHECK(NN[0]);
  CHECK_FALSE(NN[1]);
  CHECK(NN[2]);
  CHECK_FALSE(NN[3]);
  CHECK_FALSE(NN[4]);
  // ...
  CHECK(NN[23]);
  CHECK_FALSE(NN[24]);
  delete[] NN;
  NN = NULL;
  delete[] Tss__;
  Tss__ = NULL;
}
