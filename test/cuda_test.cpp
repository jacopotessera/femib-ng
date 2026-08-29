#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../src/mesh/mesh.hpp"
#include "../src/types/types.hpp"
#include "cuda.h"
#include "name_reporter.h"
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
  CHECK(femib::cuda::getStackSize() == FEMIB_CUDA_STACK_SIZE);
  // the CUDA driver rounds up to its own alignment
  CHECK(femib::cuda::getHeapSize() >= FEMIB_CUDA_HEAP_SIZE * sizeof(double));
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
  CHECK(NNN[0] == 7);
  CHECK(NNN[1] == 7);
  CHECK(NNN[2] == 7);

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

  // TODO add a better interface in cuda, that copies and then deletes and just
  //  gives the result
  delete[] NN;
  NN = NULL;
  delete[] T;
  T = NULL;

  CHECK(NNN[0] == 7);
  CHECK(NNN[1] == 7);
  CHECK(NNN[2] == 7);
}

// TODO parallel_accurate is not really accurate...
TEST_CASE("testing serial_accurate(CPU) vs parallel_accurate(GPU)") {
  mesh.init();
  femib::types::box<float, 2> boxx =
      femib::mesh::lin_spaced<float, 2>(box, delta);
  int size_T = mesh.N.size();

  bool Nser[boxx.size() * size_T];
  femib::cuda::serial_accurate<float, 2>(boxx.data(), boxx.size(),
                                         mesh.N.data(), size_T, Nser);

  femib::types::dtrian_<float, 2> *T =
      femib::types::vector_dtrian2pointer_dtrian_<float, 2>(mesh.N);
  femib::types::dtrian_<float, 2> *devT =
      femib::cuda::copyToDevice<femib::types::dtrian_<float, 2>>(T, size_T);
  femib::types::dvec<float, 2> *devX =
      femib::cuda::copyToDevice<femib::types::dvec<float, 2>>(boxx.data(),
                                                              boxx.size());
  bool Npar_init[boxx.size() * size_T];
  bool *devN = femib::cuda::copyToDevice<bool>(Npar_init, boxx.size() * size_T);
  femib::cuda::parallel_accurate<float, 2>(devX, boxx.size(), devT, size_T,
                                           devN);
  bool *Npar = femib::cuda::copyToHost<bool>(devN, boxx.size() * size_T);

  auto shared_vertices = [&](int n1, int n2) {
    int count = 0;
    for (int a = 0; a < 3; ++a) {
      for (int b = 0; b < 3; ++b) {
        if (mesh.T[n1](a) == mesh.T[n2](b)) {
          ++count;
        }
      }
    }
    return count;
  };

  int gpu_true_cpu_false = 0;
  int cpu_only_true = 0;
  int cpu_only_true_gpu_adjacent = 0;
  int cpu_only_true_gpu_elsewhere = 0;
  int cpu_only_true_gpu_lost = 0;
  for (int k = 0; k < boxx.size() * size_T; ++k) {
    if (Npar[k] && !Nser[k]) {
      ++gpu_true_cpu_false;
    }
    if (Nser[k] && !Npar[k]) {
      ++cpu_only_true;
      int i = k / size_T;
      int n = k % size_T;
      int gpu_match = -1;
      for (int n2 = 0; n2 < size_T; ++n2) {
        if (Npar[i * size_T + n2]) {
          gpu_match = n2;
          break;
        }
      }
      if (gpu_match == -1) {
        ++cpu_only_true_gpu_lost;
        std::cerr << "[cpu_only_true] point " << i << " (" << boxx[i](0) << ", "
                  << boxx[i](1) << ") cpu triangle " << n
                  << " -- gpu: unclassified in every triangle" << std::endl;
      } else {
        int shared = shared_vertices(n, gpu_match);
        if (shared >= 2) {
          ++cpu_only_true_gpu_adjacent;
        } else {
          ++cpu_only_true_gpu_elsewhere;
        }
        std::cerr << "[cpu_only_true] point " << i << " (" << boxx[i](0) << ", "
                  << boxx[i](1) << ") cpu triangle " << n << " -- gpu triangle "
                  << gpu_match << " (shared vertices " << shared << ")"
                  << std::endl;
      }
    }
  }
  CHECK(gpu_true_cpu_false == 0);
  CHECK(cpu_only_true <= 20);
  CAPTURE(cpu_only_true_gpu_adjacent);
  CHECK(cpu_only_true_gpu_lost == 0);
  CHECK(cpu_only_true_gpu_elsewhere == 0);

  delete[] Npar;
  Npar = NULL;
  delete[] T;
  T = NULL;
}