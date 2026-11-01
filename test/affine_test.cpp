#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "affine.hpp"
#include "name_reporter.h"
#include <doctest/doctest.h>

const float EPSILON = std::numeric_limits<float>::epsilon();

using namespace femib::affine;

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

// TODO use meaningful test descriptions
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
      dvec transformed_point = reference_to_physical<float, 2>(triangle, point);
      dvec inv_transf_point =
          physical_to_reference<float, 2>(triangle, transformed_point);
      dmat jacobian_matrix = jacobian<float, 2>(triangle);
      dmat inverse_jacobian_matrix = inverse_jacobian<float, 2>(triangle);

      CHECK_LE((transformed_point - linear_combination(triangle, c)).norm(),
               EPSILON * transformed_point.norm());
      CHECK_LE((inv_transf_point - point).norm(), EPSILON * 1);
      CHECK_LE((jacobian_matrix * inverse_jacobian_matrix).determinant() - 1,
               EPSILON * 1);
    }
  }
}

TEST_CASE("testing abs_jacobian_determinant") {
  femib::types::dtrian<float, 2> t = {{0.0, 0.0}, {2.0, 0.0}, {0.0, 2.0}};
  // abs_jacobian_determinant must be twice the triangle area
  float area_scaling = femib::affine::abs_jacobian_determinant<float, 2>(t);
  CHECK_EQ(area_scaling, doctest::Approx(4.0));
}

TEST_CASE("testing abs_jacobian_determinant: degenerate triangle") {
  femib::types::dtrian<float, 2> degenerate = {
      {0.0, 0.0}, {1.0, 0.0}, {2.0, 0.0}};
  float area_scaling =
      femib::affine::abs_jacobian_determinant<float, 2>(degenerate);
  CHECK_EQ(area_scaling, doctest::Approx(0.0));
}
