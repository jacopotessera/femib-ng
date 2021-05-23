#include <benchmark/benchmark.h>

#include "../src/affine/affine.hpp"
#include <Eigen/LU>
#include <algorithm>
#include <assert.h>
#include <cmath>

typedef femib::types::dvec<float, 2> dvec;
typedef femib::types::dtrian<float, 2> dtrian;
typedef femib::types::dmat<float, 2> dmat;

const std::vector<dvec> base_triangle = {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}};

dvec linear_combination(const std::vector<dvec> &triangle,
                        const std::vector<float> &a) {
  return (1 / (a[0] + a[1] + a[2])) *
         (a[0] * triangle[0] + a[1] * triangle[1] + a[2] * triangle[2]);
};

dvec base_linear_combination(const std::vector<float> &a) {
  return linear_combination(base_triangle, a);
};

const float EPSILON = std::numeric_limits<float>::epsilon();

const std::vector<std::vector<float>> coefficients = {
      {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {0.0, 0.5, 0.5},
      {0.3, 0.4, 0.3}, {1.3, 5.4, 3.3}, {0.7, 0.2, 0.1}};
const std::vector<std::vector<dvec>> triangles = {
      {{1.0, 2.0}, {3.0, 0.0}, {2.0, 4.0}},
      {{12.0, 2.0}, {33.0, 10.0}, {-22.0, 34.3}}};

const std::vector<dvec> triangle = triangles[0];
const std::vector<float> c = coefficients[0];
const dvec point = base_linear_combination(c);

const dvec transformed_point = affine<float, 2>(triangle, point);
const dvec inv_transf_point = affine_inv<float, 2>(triangle, transformed_point);
const dmat affine_matrix = affineB<float, 2>(triangle);
const dmat inv_affine_matrix = affineBinv<float, 2>(triangle);

dmat B = Eigen::Matrix<float, 2, 2>::Identity();

static void eye(benchmark::State &state) {
  for (auto _ : state) {
    B.determinant();
    }
}

BENCHMARK(eye);

static void matrix_times_inv_matrix_(benchmark::State &state) {
  for (auto _ : state) {
    (affine_matrix * affine_matrix.inverse()).determinant();
    }
}

BENCHMARK(matrix_times_inv_matrix_);

static void matrix_times_inv_matrix__(benchmark::State &state) {
  for (auto _ : state) {
    (affine_matrix.inverse() * affine_matrix).determinant();
    }
}

BENCHMARK(matrix_times_inv_matrix__);

static void matrix_times_inv_matrix(benchmark::State &state) {
  for (auto _ : state) {
    (affine_matrix * inv_affine_matrix).determinant();
    }
}

BENCHMARK(matrix_times_inv_matrix);

BENCHMARK_MAIN();
