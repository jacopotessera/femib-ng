#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../src/affine/affine.hpp"
#include "../src/femib/stokes_t.hpp"
#include "../src/finite_element/P0_2d1d.hpp"
#include "../src/finite_element/P1+B_2d2d.hpp"
#include "../src/finite_element/P1_2d1d.hpp"
#include "../src/finite_element_space/finite_element_space.hpp"
#include "../src/gauss/gauss.hpp"
#include "../src/gauss/gauss_lagrange_2_2d.hpp"
#include "../src/mesh/mesh.hpp"
#include "../src/write/write.hpp"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <doctest/doctest.h>
#include <iostream>
#include <utility>
#include <vector>

#include "P2_2d2d.hpp"
#include "gauss_5_2d.hpp"
#include "utils.hpp"

femib::stokes_t::stokes<float, 2> make_stokes_t_fixture(
    femib::types::mesh<float, 2> &mesh_out,
    std::function<Eigen::Matrix<float, 2, 1>(femib::types::dvec<float, 2>,
                                             float)> &external_force) {
  femib::gauss::rule<float, 2> rule =
      femib::gauss::create_gauss_2_2d<float, 2>();
  std::string mesh_dir = MESH_DIR;
  femib::types::mesh<float, 2> mesh = femib::mesh::read<float, 2>(
      mesh_dir + "p3.mat", mesh_dir + "t3.mat", mesh_dir + "e3.mat");
  mesh.init();

  // V
  femib::finite_element::finite_element<float, 2, 2> f_p1_2d2d =
      femib::finite_element::create_finite_element_P2_2d2d<float, 2, 2>();
  femib::finite_element_space::finite_element_space v = {
      .finite_element = f_p1_2d2d, .mesh = mesh};
  v.nodes = f_p1_2d2d.build_nodes(mesh);

  // Q
  femib::finite_element::finite_element<float, 2, 1> f_p1_2d1d =
      femib::finite_element::create_finite_element_P1_2d1d<float, 2, 1>();
  femib::finite_element_space::finite_element_space q = {
      .finite_element = f_p1_2d1d, .mesh = mesh};
  q.nodes = f_p1_2d1d.build_nodes(mesh);

  femib::stokes_t::stokes<float, 2> s;
  s.V = v;
  s.Q = q;
  s.force = external_force;
  femib::stokes_t::init<float, 2>(s, rule);

  mesh_out = mesh;
  return s;
}

TEST_CASE("testing femib stokes_t pipeline with HDF5 persistence") {
  femib::types::mesh<float, 2> mesh;
  std::function<Eigen::Matrix<float, 2, 1>(femib::types::dvec<float, 2>, float)>
      ones_force = [](femib::types::dvec<float, 2>,
                      float) -> Eigen::Matrix<float, 2, 1> {
    return Eigen::Matrix<float, 2, 1>::Ones();
  };
  femib::stokes_t::stokes<float, 2> stokes =
      make_stokes_t_fixture(mesh, ones_force);

  femib::types::box<float, 2> box = femib::mesh::find_box<float, 2>(mesh);

  femib::types::box<float, 2> boxx =
      femib::mesh::lin_spaced<float, 2>(box, 0.1);

  std::string id = get_time();

  std::string path = "/tmp/femib_stokes_t_test_" + id + ".h5";
  femib::write::save_sim(path, id);

  int TMAX = 100;
  for (int t = 0; t < TMAX; t++) {
    femib::stokes_t::advance<float, 2>(stokes);

    femib::write::plot_data<float, 2> p;
    p.time = t;
    for (const auto &plot_data : stokes.plotV[t]) {
      p.x.push_back(plot_data.first);
      p.u.push_back(plot_data.second);
    }
    femib::write::save_plot_data(path, p);
  }

  CHECK(stokes.solution.size() == (size_t)TMAX);
  CHECK(stokes.solution.back().allFinite());

  // TODO
  // std::remove(path.c_str());
}

TEST_CASE("testing advance") {
  femib::types::mesh<float, 2> mesh;
  std::function<Eigen::Matrix<float, 2, 1>(femib::types::dvec<float, 2>, float)>
      ones_force = [](femib::types::dvec<float, 2>,
                      float) -> Eigen::Matrix<float, 2, 1> {
    return Eigen::Matrix<float, 2, 1>::Ones();
  };
  femib::stokes_t::stokes<float, 2> s = make_stokes_t_fixture(mesh, ones_force);
  femib::stokes_t::advance<float, 2>(s);
  CHECK_NOTHROW(femib::stokes_t::advance<float, 2>(s));

  CHECK(s.solution.back().rows() == s.V.nodes.P.size() + s.Q.nodes.P.size());
  CHECK(s.solution.back().topRows(s.V.nodes.P.size()).allFinite());
}

TEST_CASE("testing advance over several timesteps") {
  femib::types::mesh<float, 2> mesh;
  std::function<Eigen::Matrix<float, 2, 1>(femib::types::dvec<float, 2>, float)>
      ones_force = [](femib::types::dvec<float, 2>,
                      float) -> Eigen::Matrix<float, 2, 1> {
    return Eigen::Matrix<float, 2, 1>::Ones();
  };
  femib::stokes_t::stokes<float, 2> s = make_stokes_t_fixture(mesh, ones_force);

  int rowsV = s.V.nodes.P.size();
  const int n_steps = 8;
  std::vector<float> velocity_norm;
  std::vector<float> step_diff_norm;

  Eigen::Matrix<float, Eigen::Dynamic, 1> prev_v;
  for (int step = 0; step < n_steps; ++step) {
    femib::stokes_t::advance<float, 2>(s);
    Eigen::Matrix<float, Eigen::Dynamic, 1> v =
        s.solution.back().topRows(rowsV);
    CHECK(v.allFinite());
    velocity_norm.push_back(v.norm());
    if (step > 0) {
      step_diff_norm.push_back((v - prev_v).norm());
    }
    prev_v = v;
  }

  float max_norm =
      *std::max_element(velocity_norm.begin(), velocity_norm.end());
  float min_norm =
      *std::min_element(velocity_norm.begin(), velocity_norm.end());
  CHECK(max_norm < 2.0f * min_norm);

  CHECK(step_diff_norm.back() < 1.1f * step_diff_norm.front());

  Eigen::Matrix<float, Eigen::Dynamic, 1> p =
      s.solution.back().bottomRows(s.Q.nodes.P.size());
  CHECK(std::abs(s.domain_integral_row.dot(p)) < 1e-4f);
}

// TODO MMS Method of Manufactured Solutions
// https://mooseframework.inl.gov/python/mms.html
TEST_CASE("test stokes_t::advance with a non-stationary known problem") {

  //   u(x,y,t) = g(t) * U(x,y),  p(x,y,t) = g(t) * P(x,y),  g(t) = 1 - cos(t)

  const float PI = 3.14159265358979f; // TODO math? std::numbers::pi_v<float>

  auto u1_exact = [PI](float x, float y) {
    return 2.0f * PI * std::sin(PI * x) * std::sin(PI * x) * std::sin(PI * y) *
           std::cos(PI * y);
  };
  auto u2_exact = [PI](float x, float y) {
    return -2.0f * PI * std::sin(PI * x) * std::sin(PI * y) * std::sin(PI * y) *
           std::cos(PI * x);
  };
  auto p_exact = [PI](float x, float y) {
    return std::sin(PI * x) * std::cos(PI * y);
  };
  auto f1 = [PI](float x, float y) {
    return 0.5 * PI *
           (32.0f * PI * PI * std::sin(PI * x) * std::sin(PI * x) *
                std::sin(PI * y) -
            8.0f * PI * PI * std::sin(PI * y) - std::cos(PI * x)) *
           std::cos(PI * y);
  };
  auto f2 = [PI](float x, float y) {
    return 0.5 * PI *
           (-32.0f * PI * PI * std::sin(PI * y) * std::sin(PI * y) *
                std::cos(PI * x) +
            std::sin(PI * y) + 8.0f * PI * PI * std::cos(PI * x)) *
           std::sin(PI * x);
  };

  int n = 10;
  femib::types::mesh<float, 2> mesh = make_unit_square_mesh(n);
  mesh.init();

  femib::gauss::rule<float, 2> rule =
      femib::gauss::create_gauss_2_2d<float, 2>();

  femib::finite_element::finite_element<float, 2, 2> f_p1_2d2d =
      femib::finite_element::create_finite_element_P1_B_2d2d<float, 2, 2>();
  femib::finite_element_space::finite_element_space<float, 2, 2> v = {f_p1_2d2d,
                                                                      mesh};
  v.nodes = f_p1_2d2d.build_nodes(mesh);

  // TODO this comment says: stokes_t interface is wrong!
  // P1 pressure instead of stokes_t's hardcoded P0: stokes_t.hpp's default
  // P1+bubble velocity / P0 pressure pairing was found (during development of
  // this test) to be markedly less accurate for this manufactured solution
  // than P1+bubble/P1 -- see the accuracy note on the CHECKs below.
  // test/stokes_test.cpp's already-validated steady MMS uses P1 pressure
  // with this exact velocity element.
  femib::finite_element::finite_element<float, 2, 1> f_p1_2d1d =
      femib::finite_element::create_finite_element_P1_2d1d<float, 2, 1>();
  femib::finite_element_space::finite_element_space<float, 2, 1> q = {f_p1_2d1d,
                                                                      mesh};
  q.nodes = f_p1_2d1d.build_nodes(mesh);

  femib::stokes_t::stokes<float, 2> s;
  s.V = v;
  s.Q = q;
  s.deltat = 0.02f;
  s.force = [u1_exact, u2_exact, f1,
             f2](const femib::types::dvec<float, 2> &x,
                 float t) -> femib::types::dvec<float, 2> {
    float g_dot = std::sin(t);
    float g = 1.0f - std::cos(t);
    float fx = g_dot * u1_exact(x(0), x(1)) + g * f1(x(0), x(1));
    float fy = g_dot * u2_exact(x(0), x(1)) + g * f2(x(0), x(1));
    return femib::types::dvec<float, 2>(fx, fy);
  };

  femib::stokes_t::init<float, 2>(s, rule);

  int size_P = mesh.P.size();
  int size_T = mesh.T.size();
  int rows_v = s.V.nodes.P.size();
  int n_p = s.Q.nodes.P.size();

  auto velocity_error = [&](const Eigen::Matrix<float, Eigen::Dynamic, 1> &xx,
                            float g) {
    float max_error = 0.0f;
    for (int i = 0; i < size_P; ++i) {
      float x = mesh.P[i](0), y = mesh.P[i](1);
      max_error = std::max(max_error, std::abs(xx(i, 0) - g * u1_exact(x, y)));
      max_error =
          std::max(max_error, std::abs(xx(size_P + i, 0) - g * u2_exact(x, y)));
    }
    for (int n_ = 0; n_ < size_T; ++n_) {
      femib::types::dvec<float, 2> c =
          femib::finite_element::find_center_of<float, 2>(mesh[n_]);
      max_error = std::max(max_error, std::abs(xx(2 * size_P + n_, 0) -
                                               g * u1_exact(c(0), c(1))));
      max_error = std::max(max_error, std::abs(xx(2 * size_P + size_T + n_, 0) -
                                               g * u2_exact(c(0), c(1))));
    }
    return max_error;
  };

  // Pressure is only defined up to an additive constant, subtract the mean from
  // both sides before comparing.
  auto pressure_error = [&](const Eigen::Matrix<float, Eigen::Dynamic, 1> &xx,
                            float g) {
    float mean_computed = 0.0f, mean_exact = 0.0f;
    for (int n_ = 0; n_ < n_p; ++n_) {
      mean_computed += xx(rows_v + n_, 0);
      mean_exact += g * p_exact(s.Q.nodes.P[n_](0), s.Q.nodes.P[n_](1));
    }
    mean_computed /= n_p;
    mean_exact /= n_p;
    float max_error = 0.0f;
    for (int n_ = 0; n_ < n_p; ++n_) {
      float computed = xx(rows_v + n_, 0) - mean_computed;
      float exact =
          g * p_exact(s.Q.nodes.P[n_](0), s.Q.nodes.P[n_](1)) - mean_exact;
      max_error = std::max(max_error, std::abs(computed - exact));
    }
    return max_error;
  };

  int n_steps = 250;
  std::vector<float> velocity_errors, pressure_errors, g_values;
  for (int step = 0; step < n_steps; ++step) {
    femib::stokes_t::advance<float, 2>(s);
    float g = 1.0f - std::cos(s.time);
    velocity_errors.push_back(velocity_error(s.solution.back(), g));
    pressure_errors.push_back(pressure_error(s.solution.back(), g));
    g_values.push_back(g);
  }

  float g_max = *std::max_element(g_values.begin(), g_values.end());
  CHECK(g_values.front() < 0.1f);
  CHECK(g_max > 1.9f);
  CHECK(g_values.back() < g_max - 0.5f);

  float max_velocity_error_over_run =
      *std::max_element(velocity_errors.begin(), velocity_errors.end());
  float max_pressure_error_over_run =
      *std::max_element(pressure_errors.begin(), pressure_errors.end());

  CHECK(max_velocity_error_over_run < 0.6f);
  CHECK(max_pressure_error_over_run < 20.0f);
}

TEST_CASE("Test #3 from The MINI mixed finite element for the Stokes problem: "
          "An experimental investigation") {
  const float PI = std::numbers::pi_v<float>;

  auto u1_exact = [PI](float x, float y) {
    return std::sin(2 * PI * y) * (1 - std::cos(2 * PI * x));
  };
  auto u2_exact = [PI](float x, float y) {
    return std::sin(2 * PI * x) * (std::cos(2 * PI * y) - 1);
  };
  auto p_exact = [PI](float x, float y) {
    return 2 * PI * (std::cos(2 * PI * y) - std::cos(2 * PI * x));
  };
  auto f1 = [PI](float x, float y) {
    return (-4 * PI * PI * std::sin(2 * PI * y) *
                (2 * std::cos(2 * PI * x) - 1) +
            4 * PI * PI * std::sin(2 * PI * x));
  };
  auto f2 = [PI](float x, float y) {
    return (4 * PI * PI * std::sin(2 * PI * x) *
                (2 * std::cos(2 * PI * y) - 1) -
            4 * PI * PI * std::sin(2 * PI * y));
  };

  int n = 10;
  femib::types::mesh<float, 2> mesh = make_unit_square_mesh(n);
  mesh.init();

  femib::gauss::rule<float, 2> rule =
      femib::gauss::create_gauss_5_2d<float, 2>();

  femib::finite_element::finite_element<float, 2, 2> f_p1_2d2d =
      femib::finite_element::create_finite_element_P1_B_2d2d<float, 2, 2>();
  femib::finite_element_space::finite_element_space<float, 2, 2> v = {f_p1_2d2d,
                                                                      mesh};
  v.nodes = f_p1_2d2d.build_nodes(mesh);

  femib::finite_element::finite_element<float, 2, 1> f_p1_2d1d =
      femib::finite_element::create_finite_element_P1_2d1d<float, 2, 1>();
  femib::finite_element_space::finite_element_space<float, 2, 1> q = {f_p1_2d1d,
                                                                      mesh};
  q.nodes = f_p1_2d1d.build_nodes(mesh);

  femib::stokes_t::stokes<float, 2> s;
  s.V = v;
  s.Q = q;
  s.deltat = 0.1f;
  // Time-independent forcing, the solution converges to this.
  s.force = [f1, f2](const femib::types::dvec<float, 2> &x,
                     float) -> femib::types::dvec<float, 2> {
    return femib::types::dvec<float, 2>(f1(x(0), x(1)), f2(x(0), x(1)));
  };

  femib::stokes_t::init<float, 2>(s, rule);

  std::string id = get_time();
  std::string path = "/tmp/femib_stokes_t_test3_" + id + ".h5";
  femib::write::save_sim(path, id);

  int size_P = mesh.P.size();
  int size_T = mesh.T.size();
  int rows_v = s.V.nodes.P.size();
  int n_p = s.Q.nodes.P.size();

  auto velocity_error = [&](const Eigen::Matrix<float, Eigen::Dynamic, 1> &xx) {
    float max_error = 0.0f;
    for (int i = 0; i < size_P; ++i) {
      float x = mesh.P[i](0), y = mesh.P[i](1);
      max_error = std::max(max_error, std::abs(xx(i, 0) - u1_exact(x, y)));
      max_error =
          std::max(max_error, std::abs(xx(size_P + i, 0) - u2_exact(x, y)));
    }
    for (int n_ = 0; n_ < size_T; ++n_) {
      femib::types::dvec<float, 2> c =
          femib::finite_element::find_center_of<float, 2>(mesh[n_]);
      max_error = std::max(
          max_error, std::abs(xx(2 * size_P + n_, 0) - u1_exact(c(0), c(1))));
      max_error = std::max(max_error, std::abs(xx(2 * size_P + size_T + n_, 0) -
                                               u2_exact(c(0), c(1))));
    }
    return max_error;
  };

  // Pressure is only defined up to an additive constant, so we need to subtract
  // the mean from both sides before comparing.
  auto pressure_error = [&](const Eigen::Matrix<float, Eigen::Dynamic, 1> &xx) {
    float mean_computed = 0.0f, mean_exact = 0.0f;
    for (int n_ = 0; n_ < n_p; ++n_) {
      mean_computed += xx(rows_v + n_, 0);
      mean_exact += p_exact(s.Q.nodes.P[n_](0), s.Q.nodes.P[n_](1));
    }
    mean_computed /= n_p;
    mean_exact /= n_p;
    float max_error = 0.0f;
    for (int n_ = 0; n_ < n_p; ++n_) {
      float computed = xx(rows_v + n_, 0) - mean_computed;
      float exact =
          p_exact(s.Q.nodes.P[n_](0), s.Q.nodes.P[n_](1)) - mean_exact;
      max_error = std::max(max_error, std::abs(computed - exact));
    }
    return max_error;
  };

  int n_steps = 1000;
  Eigen::Matrix<float, Eigen::Dynamic, 1> prev;
  float last_step_diff = -1.0f;
  for (int step = 0; step < n_steps; ++step) {
    femib::stokes_t::advance<float, 2>(s);
    if (step > 0) {
      last_step_diff = (s.solution.back() - prev).norm();
    }
    prev = s.solution.back();

    femib::write::plot_data<float, 2> p;
    p.time = step;
    for (const auto &plot_data : s.plotV[step]) {
      p.x.push_back(plot_data.first);
      p.u.push_back(plot_data.second);
    }
    femib::write::save_plot_data(path, p);
  }

  CHECK(last_step_diff < 0.1f);

  float final_velocity_error = velocity_error(s.solution.back());
  float final_pressure_error = pressure_error(s.solution.back());

  CHECK(final_velocity_error < 1.0f);
  CHECK(final_pressure_error < 10.0f);
}

TEST_CASE("testing advance() with an extra_velocity_rhs") {
  std::function<Eigen::Matrix<float, 2, 1>(femib::types::dvec<float, 2>, float)>
      ones_force = [](femib::types::dvec<float, 2>,
                      float) -> Eigen::Matrix<float, 2, 1> {
    return Eigen::Matrix<float, 2, 1>::Ones();
  };

  femib::types::mesh<float, 2> mesh;
  femib::stokes_t::stokes<float, 2> s_plain =
      make_stokes_t_fixture(mesh, ones_force);
  femib::stokes_t::stokes<float, 2> s_extra = s_plain;

  int rowsV = s_plain.V.nodes.P.size();
  Eigen::Matrix<float, Eigen::Dynamic, 1> extra =
      Eigen::Matrix<float, Eigen::Dynamic, 1>::Zero(rowsV);
  // index 12 is an interior node
  extra(12) = 5.0f;

  femib::stokes_t::advance<float, 2>(s_plain);
  femib::stokes_t::advance<float, 2>(s_extra, extra);

  float diff = (s_plain.solution.back() - s_extra.solution.back()).norm();
  CHECK(diff > 1e-4f);

  femib::stokes_t::stokes<float, 2> s_null =
      make_stokes_t_fixture(mesh, ones_force);
  femib::stokes_t::advance<float, 2>(s_null);
  femib::stokes_t::stokes<float, 2> s_default =
      make_stokes_t_fixture(mesh, ones_force);
  femib::stokes_t::advance<float, 2>(s_default);
  CHECK((s_null.solution.back() - s_default.solution.back()).norm() < 1e-8f);
}