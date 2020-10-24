#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "affine.hpp"
#include <algorithm>
#include <doctest/doctest.h>

std::vector<dvec> base_triangle = {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}};

dvec linear_combination(const std::vector<dvec> &triangle,
                        const std::vector<float> &a) {
  return (1 / (a[0] + a[1] + a[2])) *
         (a[0] * triangle[0] + a[1] * triangle[1] + a[2] * triangle[2]);
};

dvec base_linear_combination(const std::vector<float> &a) {
  return linear_combination(base_triangle, a);
};

float EPSILON = std::numeric_limits<float>::epsilon();

TEST_CASE("testing affine") {

  std::vector<std::vector<float>> coeff = {
      {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {0.0, 0.5, 0.5},
      {0.3, 0.4, 0.3}, {1.3, 5.4, 3.3}, {0.7, 0.2, 0.1}};
  std::vector<std::vector<dvec>> triangles = {
      {{1.0, 2.0}, {3.0, 0.0}, {2.0, 4.0}},
      {{12.0, 2.0}, {33.0, 10.0}, {-22.0, 34.3}}};

  for (std::vector<dvec> triangle : triangles) {
    for (std::vector<float> c : coeff) {
      dvec point = base_linear_combination(c);
      dvec transformed_point = affine(triangle, point);
      dvec inv_transformed_point = affine_inv(triangle, transformed_point);
      CHECK((transformed_point - linear_combination(triangle, c)).norm() <=
            EPSILON * transformed_point.norm());
      CHECK((inv_transformed_point - point).norm() <= EPSILON * 1);
    }
  }
}
