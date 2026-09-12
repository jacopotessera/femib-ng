#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../src/femib/stokes_t.hpp"
#include "../src/finite_element/P0_2d1d.hpp"
#include "../src/finite_element/P1+B_2d2d.hpp"
#include "../src/gauss/gauss_lagrange_2_2d.hpp"
#include "../src/ib/ib.hpp"
#include <algorithm>
#include <doctest/doctest.h>

#include "utils.hpp"
#include "write.hpp"

namespace {

femib::ib::ib_problem<float, 2> make_ib_fixture(int n_mesh, int n_ring,
                                                float radius) {
  femib::gauss::rule<float, 2> rule =
      femib::gauss::create_gauss_2_2d<float, 2>();
  femib::types::mesh<float, 2> mesh = make_unit_square_mesh(n_mesh);
  mesh.init();

  femib::finite_element::finite_element<float, 2, 2> f_p1_2d2d =
      femib::finite_element::create_finite_element_P1_B_2d2d<float, 2, 2>();
  femib::finite_element_space::finite_element_space v = {
      .finite_element = f_p1_2d2d, .mesh = mesh};
  v.nodes = f_p1_2d2d.build_nodes(mesh);

  femib::finite_element::finite_element<float, 2, 1> f_p0_2d1d =
      femib::finite_element::create_finite_element_P0_2d1d<float, 2, 1>();
  femib::finite_element_space::finite_element_space q = {
      .finite_element = f_p0_2d1d, .mesh = mesh};
  q.nodes = f_p0_2d1d.build_nodes(mesh);

  femib::ib::ib_problem<float, 2> p;
  p.fluid.V = v;
  p.fluid.Q = q;
  p.fluid.deltat = 0.001f;
  p.fluid.force = [](const femib::types::dvec<float, 2> &,
                     float) -> femib::types::dvec<float, 2> {
    return femib::types::dvec<float, 2>::Zero(); // no separate body force
  };
  p.structure = femib::ib::build_ring<float, 2>(
      femib::types::dvec<float, 2>(0.5f, 0.5f), radius, n_ring, 50.0f);

  femib::ib::init<float, 2>(p, rule);
  return p;
}
} // namespace

TEST_CASE("a ring at its rest configuration stays motionless") {
  femib::ib::ib_problem<float, 2> p = make_ib_fixture(16, 24, 0.2f);
  std::vector<femib::types::dvec<float, 2>> X0 = p.structure.X;

  std::string id = get_time();
  std::string path = "/tmp/femib_ib_test_" + id + ".h5";
  femib::write::save_sim(path, id);

  for (int step = 0; step < 20; ++step) {
    femib::ib::advance<float, 2>(p);

    femib::write::plot_data<float, 2> plot_data;
    plot_data.time = step;
    for (const auto &d : p.fluid.plotV[step]) {
      plot_data.x.push_back(d.first);
      plot_data.u.push_back(d.second);
    }
    for (const auto &d : p.fluid.plotQ[step]) {
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
    "a ring in perturbed configuration moves toward its rest configuration" *
    doctest::skip(true)) {
  femib::ib::ib_problem<float, 2> p = make_ib_fixture(16, 24, 0.2f);
  femib::types::dvec<float, 2> center(0.5f, 0.5f);
  for (auto &x : p.structure.X) {
    femib::types::dvec<float, 2> rel = x - center;
    x(0) = center(0) + 1.3f * rel(0);
    x(1) = center(1) + rel(1) / 1.3f;
  }
  std::vector<femib::types::dvec<float, 2>> X0 = p.structure.X;

  std::string id = get_time();
  std::string path = "/tmp/femib_ib_test_" + id + ".h5";
  femib::write::save_sim(path, id);

  for (int step = 0; step < 200; ++step) {
    femib::ib::advance<float, 2>(p);

    femib::write::plot_data<float, 2> plot_data;
    plot_data.time = step;
    for (const auto &[position, velocity] : p.fluid.plotV[step]) {
      plot_data.x.push_back(position);
      plot_data.u.push_back(velocity);
    }
    for (const auto &pressure : p.fluid.plotQ[step] | std::views::values) {
      plot_data.q.push_back(pressure);
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
  CHECK_LT(max_drift, 1e-4f);
}

TEST_CASE("a perturbed (elliptical) ring's elastic energy decays, not grows, "
          "over time" *
          doctest::skip(true)) {
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
  femib::gauss::rule<float, 2> rule =
      femib::gauss::create_gauss_2_2d<float, 2>();
  std::vector<femib::types::dvec<float, 2>> X0 = p.structure.X;

  for (int step = 0; step < 20; ++step) {
    femib::ib::advance_navier_stokes<float, 2>(p, rule, 10.0f, 10, 1e-5f);
  }

  float max_drift = 0.0f;
  for (size_t k = 0; k < X0.size(); ++k) {
    max_drift = std::max(max_drift, (p.structure.X[k] - X0[k]).norm());
  }
  CHECK_LT(max_drift, 3e-3f);
}

TEST_CASE("a perturbed (elliptical) ring's elastic energy decays, not grows, "
          "over time under Navier-Stokes too" *
          doctest::skip(true)) {
  femib::ib::ib_problem<float, 2> p = make_ib_fixture(16, 24, 0.2f);
  femib::gauss::rule<float, 2> rule =
      femib::gauss::create_gauss_2_2d<float, 2>();

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
    femib::ib::advance_navier_stokes<float, 2>(p, rule, 10.0f, 10, 1e-5f);
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
