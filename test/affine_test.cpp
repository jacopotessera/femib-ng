#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "affine.hpp"
#include <Eigen/LU>
#include <algorithm>
#include <doctest/doctest.h>

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

TEST_CASE("testing affine") {

  std::vector<std::vector<float>> coefficients = {
      {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {0.0, 0.5, 0.5},
      {0.3, 0.4, 0.3}, {1.3, 5.4, 3.3}, {0.7, 0.2, 0.1}};
  std::vector<std::vector<dvec>> triangles = {
      {{1.0, 2.0}, {3.0, 0.0}, {2.0, 4.0}},
      {{12.0, 2.0}, {33.0, 10.0}, {-22.0, 34.3}}};

  for (std::vector<dvec> triangle : triangles) {
    for (std::vector<float> c : coefficients) {
      dvec point = base_linear_combination(c);
      dvec transformed_point = affine<float, 2>(triangle, point);
      dvec inv_transf_point = affine_inv<float, 2>(triangle, transformed_point);
      dmat affine_matrix = affineB<float, 2>(triangle);
      dmat inv_affine_matrix = affineBinv<float, 2>(triangle);

      CHECK((transformed_point - linear_combination(triangle, c)).norm() <=
            EPSILON * transformed_point.norm());
      CHECK((inv_transf_point - point).norm() <= EPSILON * 1);
      CHECK((affine_matrix * inv_affine_matrix).determinant() - 1 <=
            EPSILON * 1);
    }
  }
}
