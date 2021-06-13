#include <benchmark/benchmark.h>

#include "../src/affine/affine.hpp"
#include "../src/mesh/mesh.hpp"
#include "../src/types/types.hpp"
#include "cuda.h"
#include "spdlog/spdlog.h"
#include <Eigen/LU>
#include <algorithm>
#include <assert.h>
#include <cmath>
#include <iostream>
#include <vector>

static void parallel(benchmark::State &state) {

  std::string mesh_dir = MESH_DIR;
  femib::types::mesh<float, 2> mesh = femib::mesh::read<float, 2>(
      mesh_dir + "p5.mat", mesh_dir + "t5.mat", mesh_dir + "e5.mat");
  mesh.init();
  femib::types::box<float, 2> box = femib::mesh::find_box<float, 2>(mesh);
  femib::types::box<float, 2> boxx =
      femib::mesh::lin_spaced<float, 2>(box, 0.05);
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

  spdlog::set_pattern("[%Y-%m-%dT%T] [%l] [%@@%!] %v");
  SPDLOG_INFO("[points size] found to be {}", boxx.size());
  SPDLOG_INFO("[mesh size] found to be {}", mesh.N.size());
  for (auto _ : state) {

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
    // assert(NNN[0] == 0);
  }
}

BENCHMARK(parallel);

static void serial(benchmark::State &state) {

  std::string mesh_dir = MESH_DIR;
  femib::types::mesh<float, 2> mesh = femib::mesh::read<float, 2>(
      mesh_dir + "p5.mat", mesh_dir + "t5.mat", mesh_dir + "e5.mat");
  mesh.init();
  femib::types::box<float, 2> box = femib::mesh::find_box<float, 2>(mesh);
  femib::types::box<float, 2> boxx =
      femib::mesh::lin_spaced<float, 2>(box, 0.05);

  for (auto _ : state) {

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
    // assert(NNN[0] == 0);
  }
}

BENCHMARK(serial);

BENCHMARK_MAIN();
