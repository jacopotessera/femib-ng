#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../src/mesh/mesh.hpp"
#include "../src/types/types.hpp"
#include "cuda.h"
#include "name_reporter.h"
#include "spdlog/spdlog.h"
#include <algorithm>
#include <cmath>
#include <doctest/doctest.h>
#include <iostream>
#include <memory>
#include <numbers>
#include <stdexcept>
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
  femib::cuda::set_stack_size(FEMIB_CUDA_STACK_SIZE);
  femib::cuda::set_heap_size(FEMIB_CUDA_HEAP_SIZE);
  femib::cuda::print_size();
  CHECK_EQ(femib::cuda::get_stack_size(), FEMIB_CUDA_STACK_SIZE);
  // the CUDA driver rounds up to its own alignment
  CHECK_GE(femib::cuda::get_heap_size(), FEMIB_CUDA_HEAP_SIZE * sizeof(double));
}

TEST_CASE("testing cuda copy") {
  double x = 10;
  double *X = femib::cuda::copy_to_device<double>(&x, 1);
  double *y = femib::cuda::copy_to_host<double>(X, 1);
  CHECK_EQ(x, *y);
  femib::cuda::free_device<double>(X);
  delete[] y;
}

TEST_CASE("testing cuda in_box") {
  CHECK_UNARY(femib::cuda::in_box<float, 2>(P1, T));
  CHECK_FALSE(femib::cuda::in_box<float, 2>(P2, T));
  CHECK_UNARY(femib::cuda::in_box<float, 2>(P3, T));
  CHECK_FALSE(femib::cuda::in_box<float, 2>(P4, T));
  CHECK_UNARY(femib::cuda::in_box<float, 2>(P5, T));
}

TEST_CASE("testing cuda in_triangle") {
  CHECK_UNARY(femib::cuda::in_triangle<float, 2>(P1, T));
  CHECK_FALSE(femib::cuda::in_triangle<float, 2>(P2, T));
  CHECK_UNARY(femib::cuda::in_triangle<float, 2>(P3, T));
  CHECK_FALSE(femib::cuda::in_triangle<float, 2>(P4, T));
  CHECK_FALSE(femib::cuda::in_triangle<float, 2>(P5, T));
}

TEST_CASE("testing cuda accurate") {
  CHECK_UNARY(femib::cuda::accurate<float, 2>(P1, T));
  CHECK_FALSE(femib::cuda::accurate<float, 2>(P2, T));
  CHECK_UNARY(femib::cuda::accurate<float, 2>(P3, T));
  CHECK_FALSE(femib::cuda::accurate<float, 2>(P4, T));
  CHECK_FALSE(femib::cuda::accurate<float, 2>(P5, T));
}

TEST_CASE("testing cuda serial_accurate") {
  std::vector<femib::types::dvec<float, 2>> grid =
      femib::mesh::build_uniform_grid<float, 2>(box, delta);

  bool N[grid.size() * mesh.N.size()];

  femib::cuda::serial_accurate<float, 2>(grid.data(), grid.size(),
                                         mesh.N.data(), mesh.N.size(), N);

  std::vector<int> NNN;

  for (size_t i = 0; i < grid.size(); ++i) {
    for (int n = 0; n < std::ssize(mesh.N); ++n) {
      if (N[i * mesh.N.size() + n]) {
        NNN.push_back(n);
        break;
      }
    }
  }
  CHECK_EQ(NNN[0], 7);
  CHECK_EQ(NNN[1], 7);
  CHECK_EQ(NNN[2], 7);

  for (size_t i = 0; i < grid.size(); ++i) {
    CHECK_GE(NNN[i], 0);
    CHECK_LT(NNN[i], mesh.N.size());
  }
}

TEST_CASE("testing cuda parallel_accurate") {
  std::vector<femib::types::dvec<float, 2>> grid =
      femib::mesh::build_uniform_grid<float, 2>(box, delta);

  bool N[grid.size() * mesh.N.size()];

  femib::types::dtrian_<float, 2> *T =
      femib::types::vector_dtrian2pointer_dtrian_<float, 2>(mesh.N);

  femib::types::dtrian_<float, 2> *devT =
      femib::cuda::copy_to_device<femib::types::dtrian_<float, 2>>(
          T, mesh.N.size());
  femib::types::dvec<float, 2> *devX =
      femib::cuda::copy_to_device<femib::types::dvec<float, 2>>(grid.data(),
                                                                grid.size());
  bool *devN =
      femib::cuda::copy_to_device<bool>(N, grid.size() * mesh.N.size());
  femib::cuda::parallel_accurate<float, 2>(devX, grid.size(), devT,
                                           mesh.N.size(), devN);
  bool *NN;
  NN = femib::cuda::copy_to_host<bool>(devN, grid.size() * mesh.N.size());

  std::vector<int> NNN;

  for (size_t i = 0; i < grid.size(); ++i) {
    for (int n = 0; n < std::ssize(mesh.N); ++n) {
      if (NN[i * mesh.N.size() + n]) {
        NNN.push_back(n);
        break;
      }
    }
  }

  // TODO add a better interface in cuda, that copies and then deletes and just
  //  gives the result
  femib::cuda::free_device<femib::types::dtrian_<float, 2>>(devT);
  femib::cuda::free_device<femib::types::dvec<float, 2>>(devX);
  femib::cuda::free_device<bool>(devN);
  delete[] NN;
  NN = NULL;
  delete[] T;
  T = NULL;

  CHECK_EQ(NNN[0], 7);
  CHECK_EQ(NNN[1], 7);
  CHECK_EQ(NNN[2], 7);
}

TEST_CASE("parallel_accurate checks the grid limits before launching") {
  // every GPU since compute capability 3.0: grid.x <= 2^31 - 1, grid.y <= 65535
  const size_t too_many_points = size_t{1} << 31;
  const size_t too_many_triangles = size_t{65535} * 256 + 1;
  CHECK_THROWS_AS((femib::cuda::parallel_accurate<float, 2>(
                      nullptr, too_many_points, nullptr, 1, nullptr)),
                  std::invalid_argument);
  CHECK_THROWS_AS((femib::cuda::parallel_accurate<float, 2>(
                      nullptr, 1, nullptr, too_many_triangles, nullptr)),
                  std::invalid_argument);
  CHECK_NOTHROW((femib::cuda::parallel_accurate<float, 2>(nullptr, 0, nullptr,
                                                          10, nullptr)));
}

TEST_CASE("testing mesh find_points") {
  std::vector<femib::types::dvec<float, 2>> grid =
      femib::mesh::build_uniform_grid<float, 2>(box, delta);

  std::vector<int> NNN = femib::mesh::find_points<float, 2>(mesh, grid);

  CHECK_EQ(NNN[0], 7);
  CHECK_EQ(NNN[1], 7);
  CHECK_EQ(NNN[2], 7);
}

TEST_CASE("testing serial_accurate(CPU) vs parallel_accurate(GPU)") {
  std::vector<femib::types::dvec<float, 2>> grid =
      femib::mesh::build_uniform_grid<float, 2>(box, delta);
  int size_T = mesh.N.size();

  bool Nser[grid.size() * size_T];
  femib::cuda::serial_accurate<float, 2>(grid.data(), grid.size(),
                                         mesh.N.data(), size_T, Nser);

  femib::types::dtrian_<float, 2> *T =
      femib::types::vector_dtrian2pointer_dtrian_<float, 2>(mesh.N);
  femib::types::dtrian_<float, 2> *devT =
      femib::cuda::copy_to_device<femib::types::dtrian_<float, 2>>(T, size_T);
  femib::types::dvec<float, 2> *devX =
      femib::cuda::copy_to_device<femib::types::dvec<float, 2>>(grid.data(),
                                                                grid.size());
  bool Npar_init[grid.size() * size_T];
  bool *devN =
      femib::cuda::copy_to_device<bool>(Npar_init, grid.size() * size_T);
  femib::cuda::parallel_accurate<float, 2>(devX, grid.size(), devT, size_T,
                                           devN);
  bool *Npar = femib::cuda::copy_to_host<bool>(devN, grid.size() * size_T);

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
  for (int k = 0; k < std::ssize(grid) * size_T; ++k) {
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
        std::cerr << "[cpu_only_true] point " << i << " (" << grid[i](0) << ", "
                  << grid[i](1) << ") cpu triangle " << n
                  << " -- gpu: unclassified in every triangle" << std::endl;
      } else {
        int shared = shared_vertices(n, gpu_match);
        if (shared >= 2) {
          ++cpu_only_true_gpu_adjacent;
        } else {
          ++cpu_only_true_gpu_elsewhere;
        }
        std::cerr << "[cpu_only_true] point " << i << " (" << grid[i](0) << ", "
                  << grid[i](1) << ") cpu triangle " << n << " -- gpu triangle "
                  << gpu_match << " (shared vertices " << shared << ")"
                  << std::endl;
      }
    }
  }
  CHECK_EQ(gpu_true_cpu_false, 0);
  CHECK_EQ(cpu_only_true, 0);
  CAPTURE(cpu_only_true_gpu_adjacent);
  CHECK_EQ(cpu_only_true_gpu_lost, 0);
  CHECK_EQ(cpu_only_true_gpu_elsewhere, 0);

  femib::cuda::free_device<femib::types::dtrian_<float, 2>>(devT);
  femib::cuda::free_device<femib::types::dvec<float, 2>>(devX);
  femib::cuda::free_device<bool>(devN);
  delete[] Npar;
  Npar = NULL;
  delete[] T;
  T = NULL;
}
TEST_CASE("points lying exactly on a mesh's boundary are located") {
  using vec = femib::types::dvec<double, 2>;
  // [0,2]x[0,1] in 40x20 squares split in two; then the same, turned by 30
  // degrees so that its boundary edges aren't axis-aligned
  const int nx = 40, ny = 20;
  auto idx = [](int i, int j) { return i * (ny + 1) + j; };
  std::vector<femib::types::ditrian<2>> triangles;
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      triangles.emplace_back(idx(i, j), idx(i + 1, j), idx(i + 1, j + 1));
      triangles.emplace_back(idx(i, j), idx(i + 1, j + 1), idx(i, j + 1));
    }
  }
  const double angle = std::numbers::pi / 6;
  for (double turn : {0.0, angle}) {
    auto place = [turn](double x, double y) {
      return vec(std::cos(turn) * x - std::sin(turn) * y,
                 std::sin(turn) * x + std::cos(turn) * y);
    };
    std::vector<vec> nodes;
    for (int i = 0; i <= nx; ++i) {
      for (int j = 0; j <= ny; ++j) {
        nodes.push_back(place(2.0 * i / nx, 1.0 * j / ny));
      }
    }
    femib::types::mesh<double, 2> rectangle(nodes, triangles);

    const std::vector<vec> corners = {place(0, 0), place(2, 0), place(2, 1),
                                      place(0, 1)};
    std::vector<vec> points;
    const int n = 250;
    for (int c = 0; c < 4; ++c) {
      for (int k = 0; k < n; ++k) {
        double s = static_cast<double>(k) / n;
        points.push_back((1 - s) * corners[c] + s * corners[(c + 1) % 4]);
      }
    }
    CAPTURE(turn);

    std::vector<int> located =
        femib::mesh::find_points<double, 2>(rectangle, points);
    CHECK_EQ(std::ranges::count(located, -1), 0);

    const size_t size_T = rectangle.N.size();
    auto N = std::make_unique<bool[]>(points.size() * size_T);
    femib::cuda::serial_accurate<double, 2>(
        points.data(), points.size(), rectangle.N.data(), size_T, N.get());
    int cpu_missed = 0;
    for (size_t i = 0; i < points.size(); ++i) {
      cpu_missed +=
          std::none_of(N.get() + i * size_T, N.get() + (i + 1) * size_T,
                       [](bool found) { return found; });
    }
    CHECK_EQ(cpu_missed, 0);
  }
}
