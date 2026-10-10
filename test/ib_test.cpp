#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../src/finite_element/P0_2d1d.hpp"
#include "../src/finite_element/P1+B_2d2d.hpp"
#include "../src/gauss/gauss_2_2d.hpp"
#include "../src/ib/ib.hpp"
#include <algorithm>
#include <cmath>
#include <doctest/doctest.h>
#include <memory>
#include <vector>

#include "utils.hpp"
#include "write.hpp"

namespace {

femib::ib::ib_problem<float, 2>
make_ib_fixture(int n_mesh, int n_ring, float radius, float k_spring = 2.618f) {
  femib::gauss::rule<float, 2> rule =
      femib::gauss::create_gauss_2_2d<float, 2>();
  femib::types::mesh<float, 2> mesh = make_unit_square_mesh(n_mesh);

  femib::finite_element::finite_element<float, 2, 2> f_p1_2d2d =
      femib::finite_element::create_finite_element_P1_B_2d2d<float, 2, 2>();
  femib::finite_element_space::finite_element_space v(f_p1_2d2d, mesh);

  femib::finite_element::finite_element<float, 2, 1> f_p0_2d1d =
      femib::finite_element::create_finite_element_P0_2d1d<float, 2, 1>();
  femib::finite_element_space::finite_element_space q(f_p0_2d1d, mesh);

  femib::ib::ib_problem<float, 2> p(
      v, q, rule, /*rho=*/1.0f, /*mu=*/1.0f, /*deltat=*/0.001f,
      [](const femib::types::dvec<float, 2> &,
         float) -> femib::types::dvec<float, 2> {
        return femib::types::dvec<float, 2>::Zero(); // no separate body force
      },
      std::make_unique<femib::util::picard_solver<float>>(10, 1e-5f));
  p.structure = femib::ib::build_ring<float, 2>(
      femib::types::dvec<float, 2>(0.5f, 0.5f), radius, n_ring, k_spring);
  return p;
}

// sqrt(var(x) / var(y)) of the points: a/b for an axis-aligned ellipse
float moment_aspect(const std::vector<femib::types::dvec<float, 2>> &X) {
  femib::types::dvec<float, 2> c = femib::types::dvec<float, 2>::Zero();
  for (const auto &x : X) {
    c += x;
  }
  c /= static_cast<float>(X.size());
  float var_x = 0.0f;
  float var_y = 0.0f;
  for (const auto &x : X) {
    var_x += (x(0) - c(0)) * (x(0) - c(0));
    var_y += (x(1) - c(1)) * (x(1) - c(1));
  }
  return std::sqrt(var_x / var_y);
}

femib::types::dvec<float, 2>
extent(const std::vector<femib::types::dvec<float, 2>> &X) {
  femib::types::dvec<float, 2> lo = X[0];
  femib::types::dvec<float, 2> hi = X[0];
  for (const auto &x : X) {
    lo = lo.cwiseMin(x);
    hi = hi.cwiseMax(x);
  }
  return hi - lo;
}

float polygon_area(const std::vector<femib::types::dvec<float, 2>> &X) {
  float A = 0.0f;
  for (size_t k = 0; k < X.size(); ++k) {
    const auto &a = X[k];
    const auto &b = X[(k + 1) % X.size()];
    A += a(0) * b(1) - b(0) * a(1);
  }
  return std::abs(A) * 0.5f;
}
} // namespace

TEST_CASE("a ring at its rest configuration stays motionless") {
  femib::ib::ib_problem<float, 2> p = make_ib_fixture(16, 24, 0.2f);
  std::vector<femib::types::dvec<float, 2>> X0 = p.structure.X;

  std::string id = get_time();
  std::string path = "/tmp/femib_ib_test_" + id + ".h5";
  femib::write::save_sim(path, id, "ring", femib::write::mode::overwrite);

  for (int step = 0; step < 20; ++step) {
    femib::ib::advance<float, 2>(p);

    femib::write::plot_data<float, 2> plot_data;
    plot_data.time = step;
    for (const auto &d : p.fluid.plot_velocity(0.01f)) {
      plot_data.x.push_back(d.first);
      plot_data.u.push_back(d.second);
    }
    for (const auto &d : p.fluid.plot_pressure(0.01f)) {
      plot_data.q.push_back(d.second);
    }
    for (const auto &val : p.structure.X) {
      plot_data.X.push_back(val);
    }
    femib::write::save_plot_data(path, plot_data);
  }

  float max_drift = 0.0f;
  for (size_t k = 0; k < X0.size(); ++k) {
    max_drift = std::max(max_drift, (p.structure.X[k] - X0[k]).norm());
  }
  CHECK_LT(max_drift, 3e-3f);
}

TEST_CASE(
    "a ring in perturbed configuration moves toward its rest configuration") {
  // n_ring keeps the ring's point spacing under half the mesh spacing, so the
  // boundary doesn't leak; k_spring=100 is unstable at this mesh and deltat
  femib::ib::ib_problem<float, 2> p = make_ib_fixture(24, 96, 0.2f, 50.0f);
  femib::types::dvec<float, 2> center(0.5f, 0.5f);
  for (auto &x : p.structure.X) {
    femib::types::dvec<float, 2> rel = x - center;
    x(0) = center(0) + 1.3f * rel(0);
    x(1) = center(1) + rel(1) / 1.3f;
  }
  std::vector<femib::types::dvec<float, 2>> X0 = p.structure.X;
  const femib::types::dvec<float, 2> extent0 = extent(X0);
  const float area0 = polygon_area(X0);

  float aspect = moment_aspect(X0);
  for (int step = 1; step <= 200; ++step) {
    femib::ib::advance<float, 2>(p);
    if (step % 50 == 0) {
      const float new_aspect = moment_aspect(p.structure.X);
      CHECK_LT(new_aspect, aspect);
      aspect = new_aspect;
    }
  }

  float max_drift = 0.0f;
  femib::types::dvec<float, 2> centroid = femib::types::dvec<float, 2>::Zero();
  for (size_t k = 0; k < X0.size(); ++k) {
    max_drift = std::max(max_drift, (p.structure.X[k] - X0[k]).norm());
    centroid += p.structure.X[k];
  }
  centroid /= static_cast<float>(X0.size());
  const femib::types::dvec<float, 2> extent1 = extent(p.structure.X);

  CHECK_GT(max_drift, 0.015f);
  CHECK_LT(aspect, 1.62f); // starts at 1.3^2 = 1.69
  // shrinking alone would narrow both axes: the short one grows only if the
  // ring is rounding out
  CHECK_LT(extent1(0), 0.97f * extent0(0));
  CHECK_GT(extent1(1), 1.02f * extent0(1));
  CHECK_GT(polygon_area(p.structure.X), 0.9f * area0);
  CHECK_LT((centroid - center).norm(), 2e-3f);
}

TEST_CASE("a perturbed (elliptical) ring's elastic energy decays, not grows, "
          "over time") {
  femib::ib::ib_problem<float, 2> p = make_ib_fixture(16, 24, 0.2f);

  femib::types::dvec<float, 2> center(0.5f, 0.5f);
  for (auto &x : p.structure.X) {
    femib::types::dvec<float, 2> rel = x - center;
    x(0) = center(0) + 1.3f * rel(0);
    x(1) = center(1) + rel(1) / 1.3f;
  }

  auto E0 = femib::ib::elastic_energy<float, 2>(p.structure);
  REQUIRE_GT(E0,
             1e-6f); // sanity check: the perturbation actually stores energy

  std::vector<float> energies;
  int n_steps = 200;
  for (int step = 0; step < n_steps; ++step) {
    femib::ib::advance<float, 2>(p);
    energies.push_back(femib::ib::elastic_energy<float, 2>(p.structure));
  }

  float E_early = energies[9]; // after a short transient
  float E_late = energies[n_steps - 1];
  float E_max = *std::ranges::max_element(energies);

  CHECK_LT(E_late, E_early);
  CHECK_LT(E_max, 2.0f * E0);
  CHECK_LT(E_late, 0.985f * E0);
}

TEST_CASE("a ring at its rest configuration stays motionless under "
          "Navier-Stokes too") {
  femib::ib::ib_problem<float, 2> p = make_ib_fixture(16, 24, 0.2f);
  std::vector<femib::types::dvec<float, 2>> X0 = p.structure.X;

  for (int step = 0; step < 20; ++step) {
    femib::ib::advance_navier_stokes<float, 2>(p);
  }

  float max_drift = 0.0f;
  for (size_t k = 0; k < X0.size(); ++k) {
    max_drift = std::max(max_drift, (p.structure.X[k] - X0[k]).norm());
  }
  CHECK_LT(max_drift, 3e-3f);
}

TEST_CASE("a perturbed (elliptical) ring's elastic energy decays, not grows, "
          "over time under Navier-Stokes too") {
  femib::ib::ib_problem<float, 2> p = make_ib_fixture(16, 24, 0.2f);

  femib::types::dvec<float, 2> center(0.5f, 0.5f);
  for (auto &x : p.structure.X) {
    femib::types::dvec<float, 2> rel = x - center;
    x(0) = center(0) + 1.3f * rel(0);
    x(1) = center(1) + rel(1) / 1.3f;
  }

  auto E0 = femib::ib::elastic_energy<float, 2>(p.structure);
  REQUIRE_GT(E0, 1e-6f); // sanity: the perturbation actually stored energy

  float E_max = E0;
  int n_steps = 20;
  for (int step = 0; step < n_steps; ++step) {
    femib::ib::advance_navier_stokes<float, 2>(p);
    E_max = std::max(E_max, femib::ib::elastic_energy<float, 2>(p.structure));
  }

  CHECK_LT(E_max, 2.0f * E0);
}

TEST_CASE("advance() throws if structure.dS is left at its default") {
  femib::ib::ib_problem<float, 2> p = make_ib_fixture(16, 24, 0.2f);
  femib::ib::ring<float, 2> bad_structure;
  bad_structure.X = p.structure.X;
  bad_structure.k = p.structure.k;
  // bad_structure.dS left at its default of 0
  p.structure = bad_structure;

  CHECK_THROWS_AS((femib::ib::advance<float, 2>(p)), std::invalid_argument);
}
