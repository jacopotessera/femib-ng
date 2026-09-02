#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../src/ib/structure.hpp"
#include "name_reporter.h"
#include <cmath>
#include <doctest/doctest.h>
#include <numbers>

#include "finite_element.hpp"

TEST_CASE("make_ring places n_points evenly on a circle of the given radius") {
  femib::types::dvec<float, 2> center(0.5f, 0.5f);
  float radius = 0.2f;
  int n = 16;
  femib::ib::ring<float, 2> r =
      femib::ib::build_ring<float, 2>(center, radius, n, 1.0f);

  REQUIRE_EQ(r.X.size(), static_cast<size_t>(n));

  for (int k = 0; k < n; ++k) {
    float dist = (r.X[k] - center).norm();
    CHECK_LT(std::abs(dist - radius), 1e-5f);
  }

  float expected_chord =
      2.0f * radius *
      std::sin(std::numbers::pi_v<float> / static_cast<float>(n));
  for (int i = 0; i < n; ++i) {
    int i_next = femib::ib::get_next_point(i, n);
    float len = (r.X[i_next] - r.X[i]).norm();
    CHECK_LT(std::abs(len - expected_chord), 1e-4f);
  }

  float expected_dS =
      (2.0f * std::numbers::pi_v<float>) / static_cast<float>(n);
  CHECK_LT(std::abs(r.dS - expected_dS), 1e-4f);
}

TEST_CASE("elastic_force goes to zero as radius goes to zero") {
  femib::types::dvec<float, 2> center(0.5f, 0.5f);
  int n_points = 16;
  float k = 3.0f;
  femib::ib::ring<float, 2> r1 =
      femib::ib::build_ring<float, 2>(center, 0.2f, n_points, k);
  femib::ib::ring<float, 2> r2 =
      femib::ib::build_ring<float, 2>(center, 0.1f, n_points, k);
  femib::ib::ring<float, 2> r3 =
      femib::ib::build_ring<float, 2>(center, 0.05f, n_points, k);

  std::vector<femib::types::dvec<float, 2>> F1 =
      femib::ib::elastic_force<float, 2>(r1);

  std::vector<femib::types::dvec<float, 2>> F2 =
      femib::ib::elastic_force<float, 2>(r2);

  std::vector<femib::types::dvec<float, 2>> F3 =
      femib::ib::elastic_force<float, 2>(r3);

  for (int i = 0; i < n_points; ++i) {
    // in a circle, the force always points towards the center
    CHECK_LT((F1[i].normalized() + (r1.X[i] - center).normalized()).norm(),
             1e-5f);
    CHECK_LT((F2[i].normalized() + (r2.X[i] - center).normalized()).norm(),
             1e-5f);
    CHECK_LT((F3[i].normalized() + (r3.X[i] - center).normalized()).norm(),
             1e-5f);

    CHECK_GT(F1[i].norm(), F2[i].norm());
    CHECK_GT(F2[i].norm(), F3[i].norm());
  }
}

TEST_CASE(
    "elastic_force matches the finite-difference gradient of elastic_energy") {
  femib::types::dvec<float, 2> center(0.5f, 0.5f);
  femib::ib::ring<float, 2> r =
      femib::ib::build_ring<float, 2>(center, 0.2f, 8, 3.0f);

  // Perturb away from the rest configuration so the force is non-trivial.
  for (size_t k = 0; k < r.X.size(); ++k) {
    r.X[k](0) += 0.01f * std::sin(3.0f * static_cast<float>(k));
    r.X[k](1) += 0.01f * std::cos(2.0f * static_cast<float>(k));
  }

  std::vector<femib::types::dvec<float, 2>> F =
      femib::ib::elastic_force<float, 2>(r);

  for (size_t k = 0; k < r.X.size(); ++k) {
    for (int comp = 0; comp < 2; ++comp) {
      float h = 1e-3f;
      femib::ib::ring<float, 2> r_plus = r;
      femib::ib::ring<float, 2> r_minus = r;
      r_plus.X[k](comp) += h;
      r_minus.X[k](comp) -= h;
      float dE_dXk = (femib::ib::elastic_energy<float, 2>(r_plus) -
                      femib::ib::elastic_energy<float, 2>(r_minus)) /
                     (2 * h);
      // force = -grad(energy).
      CHECK_LT(std::abs(F[k](comp) - (-dE_dXk)), 1e-2f);
    }
  }
}
