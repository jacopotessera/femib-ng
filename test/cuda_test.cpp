#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../src/mesh/mesh.hpp"
#include "../src/types/types.hpp"
#include "cuda.h"
#include "spdlog/spdlog.h"
#include <doctest/doctest.h>
#include <iostream>
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

std::string mesh_dir = MESH_DIR;
femib::types::mesh<float, 2> mesh = femib::mesh::read<float, 2>(
    mesh_dir + "p3.mat", mesh_dir + "t3.mat", mesh_dir + "e3.mat");
femib::types::box<float, 2> box = femib::mesh::find_box<float, 2>(mesh);
float delta = 0.019;

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
  mesh.init();
  femib::types::box<float, 2> boxx =
      femib::mesh::lin_spaced<float, 2>(box, delta);

  bool N[boxx.size() * mesh.N.size()];

  femib::cuda::serial_accurate<float, 2>(boxx.data(), boxx.size(),
                                         mesh.N.data(), mesh.N.size(), N);

  std::vector<int> NNN;

  for (int i = 0; i < boxx.size(); ++i) {
    for (int n = 0; n < mesh.N.size(); ++n) {
      if (N[i * mesh.N.size() + n]) {
        NNN.push_back(n);
        break;
      }
    }
  }
  // CHECK(NNN[0] == 0);
  // CHECK(NNN[1] == 3);
  // CHECK(NNN[2] == 3);

  for (int i = 0; i < boxx.size(); ++i) {
    CHECK(NNN[i] >= 0);
    CHECK(NNN[i] < mesh.N.size());
  }
}

TEST_CASE("testing cuda parallel_accurate") {
  mesh.init();
  femib::types::box<float, 2> boxx =
      femib::mesh::lin_spaced<float, 2>(box, delta);

  bool N[boxx.size() * mesh.N.size()];

  femib::types::dtrian_<float, 2> *T =
      femib::types::vector_dtrian2pointer_dtrian_<float, 2>(mesh.N);

  femib::types::dtrian_<float, 2> *devT =
      femib::cuda::copyToDevice<femib::types::dtrian_<float, 2>>(T,
                                                                 mesh.N.size());
  femib::types::dvec<float, 2> *devX =
      femib::cuda::copyToDevice<femib::types::dvec<float, 2>>(boxx.data(),
                                                              boxx.size());
  bool *devN = femib::cuda::copyToDevice<bool>(N, boxx.size() * mesh.N.size());
  femib::cuda::parallel_accurate<float, 2>(devX, boxx.size(), devT,
                                           mesh.N.size(), devN);
  bool *NN;
  NN = femib::cuda::copyToHost<bool>(devN, boxx.size() * mesh.N.size());

  std::vector<int> NNN;

  for (int i = 0; i < boxx.size(); ++i) {
    for (int n = 0; n < mesh.N.size(); ++n) {
      if (NN[i * mesh.N.size() + n]) {
        NNN.push_back(n);
        break;
      }
    }
  }
  // CHECK(NN[0]);
  delete[] NN;
  NN = NULL;
  delete[] T;
  T = NULL;
  for (int i = 0; i < boxx.size(); ++i) {
    // SPDLOG_INFO("[({},{})] found to be in triangle {}", boxx[i](0),
    // boxx[i](1),
    //           NNN[i]);
  }
  // CHECK(NNN[0] == 0);
  // CHECK(NNN[1] == 3);
  // CHECK(NNN[2] == 3);

  for (int i = 0; i < boxx.size(); ++i) {
    CHECK(NNN[i] >= 0);
    CHECK(NNN[i] < mesh.N.size());

    std::cerr << boxx[i](0) << "\t" << boxx[i](1) << "\t" << NNN[i]
              << std::endl;
  }
}
