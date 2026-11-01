#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../femib/femib.hpp"
#include "P1_2d1d.hpp"
#include "P1_2d2d.hpp"
#include "finite_element_space.hpp"
#include "name_reporter.h"
#include <Eigen/Sparse>
#include <doctest/doctest.h>
#include <iostream>

#include "gauss_2_2d.hpp"

int global_index(const femib::types::nodes<float, 2> &nodes, int i, int n) {
  return nodes.T[n][i];
}

TEST_CASE("testing finite_element_space") {

  std::vector<Eigen::Triplet<float>> M;
  std::vector<Eigen::Triplet<float>> F;

  femib::gauss::rule<float, 2> rule =
      femib::gauss::create_gauss_2_2d<float, 2>();

  std::string mesh_dir = MESH_DIR;
  femib::types::mesh<float, 2> mesh = femib::mesh::read<float, 2>(
      mesh_dir + "p5.mat", mesh_dir + "t5.mat", mesh_dir + "e5.mat");

  femib::finite_element::finite_element<float, 2, 1> f =
      femib::finite_element::create_finite_element_P1_2d1d<float, 2, 1>();

  femib::finite_element_space::finite_element_space<float, 2, 1> s(f, mesh);

  for (int n = 0; n < std::ssize(s.mesh.T); ++n) {
    for (int i = 0; i < std::ssize(s.finite_element.base_functions); ++i) {
      for (int j = 0; j < std::ssize(s.finite_element.base_functions); ++j) {
        femib::types::dtrian<float, 2> t = s.mesh[n];
        femib::types::dmat<float, 2> Binv = femib::affine::affineBinv(t);
        femib::types::dvec<float, 2> bb = femib::affine::affineb(t);
        femib::types::F<float, 2, 1> a =
            femib::util::reference_to_physical_function<float, 2, 1>(s, i, Binv,
                                                                     bb);
        femib::types::F<float, 2, 1> b =
            femib::util::reference_to_physical_function<float, 2, 1>(s, j, Binv,
                                                                     bb);

        float m = femib::mesh::integrate<float, 2>(
            rule,
            [&a, &b](femib::types::dvec<float, 2> x) {
              femib::types::dvec<float, 2> av = a.dx(x);
              femib::types::dvec<float, 2> bv = b.dx(x);
              return av[0] * bv[0] + av[1] * bv[1];
            },
            t);
        float f_ = femib::mesh::integrate<float, 2>(
            rule,
            [&t, &f, i](femib::types::dvec<float, 2> x) {
              float a_0 =
                  -40 * x(0) * x(1) *
                  (f.base_functions[i].x(femib::affine::affineBinv(t) *
                                         (x - femib::affine::affineb(t))))[0];
              return a_0;
            },
            t);
        M.emplace_back(global_index(s.nodes, i, n), global_index(s.nodes, j, n),
                       m);
        F.emplace_back(global_index(s.nodes, i, n), 0, f_);
      }
    }
  }
  std::vector<Eigen::Triplet<float>> B;

  std::function<float(femib::types::dvec<float, 2>)> b =
      [](const femib::types::dvec<float, 2> &x) { return x(0) + x(1); };
  //[](const femib::types::dvec<float, 2> &x) { return 1.0; };
  for (int e : s.nodes.E) {
    B.emplace_back(e, 0, b(s.nodes.P[e]));
  }

  std::vector<int> not_edge_nodes;
  for (int i = 0; i < std::ssize(s.nodes.P); i++) {
    if (std::find(s.nodes.E.begin(), s.nodes.E.end(), i) == s.nodes.E.end()) {
      not_edge_nodes.push_back(i);
    }
  }

  Eigen::SparseMatrix<float> sB =
      Eigen::SparseMatrix<float>(s.nodes.P.size(), 1);
  sB.setFromTriplets(B.begin(), B.end());
  Eigen::Matrix<float, Eigen::Dynamic, 1> dB =
      Eigen::Matrix<float, Eigen::Dynamic, 1>(sB);
  Eigen::SparseMatrix<float> sM =
      Eigen::SparseMatrix<float>(s.nodes.P.size(), s.nodes.P.size());
  sM.setFromTriplets(M.begin(), M.end());
  Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> dM =
      Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>(sM);
  Eigen::Matrix<float, Eigen::Dynamic, 1> ss = dM * dB;
  Eigen::Matrix<float, Eigen::Dynamic, 1> mf;
  mf.resize(s.nodes.P.size(), 1);
  for (int i = 0; i < std::ssize(s.nodes.P); i++) {
    auto k = std::find(not_edge_nodes.begin(), not_edge_nodes.end(), i);
    mf(i, 0) = 0.0;
    if (k != not_edge_nodes.end()) {
      mf(i, 0) = ss(k - not_edge_nodes.begin(), 0);
    }
  }

  Eigen::SparseMatrix<float> sF =
      Eigen::SparseMatrix<float>(s.nodes.P.size(), 1);
  sF.setFromTriplets(F.begin(), F.end());
  Eigen::Matrix<float, Eigen::Dynamic, 1> dF =
      Eigen::Matrix<float, Eigen::Dynamic, 1>(sF);
  Eigen::Matrix<float, Eigen::Dynamic, 1> bbb = (dF - ss)(not_edge_nodes, 0);
  Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> AAA =
      (dM)(not_edge_nodes, not_edge_nodes);
  Eigen::Matrix<float, Eigen::Dynamic, 1> xxx =
      AAA.colPivHouseholderQr().solve(bbb);
  Eigen::Matrix<float, Eigen::Dynamic, 1> xx;
  xx.resize(s.nodes.P.size(), 1);

  for (int i = 0; i < std::ssize(s.nodes.P); i++) {
    xx(i, 0) = 0.0;
    auto k = std::find(not_edge_nodes.begin(), not_edge_nodes.end(), i);
    if (k != not_edge_nodes.end()) {
      xx(i, 0) = xxx(k - not_edge_nodes.begin(), 0);
    }
    auto kk = std::find(s.nodes.E.begin(), s.nodes.E.end(), i);
    if (kk != s.nodes.E.end()) {
      xx(i, 0) = dB(i, 0);
    }
  }

  CHECK_EQ(xx.size(), s.nodes.P.size());
  for (int e : s.nodes.E) {
    CHECK_EQ(xx(e, 0), doctest::Approx(b(s.nodes.P[e])));
  }
}

TEST_CASE("plot() does not read out of bounds on a non-rectangular mesh") {
  // An L-shaped (non-convex, non-rectangular) mesh, so lin_spaced's
  // rectangular bounding box necessarily includes sample points outside
  // the triangulated region.
  //
  // plot() is only ever instantiated with e == d in production code (see
  // femib::stokes<T,d>::V, a finite_element_space<T,d,d>), since its
  // interior arithmetic hard-codes F<T,d,d> and a d-dimensional result
  // vector. P1_2d1d (e=1) does not type-check through plot(), so this
  // fixture uses the matching vector element P1_2d2d (e=d=2) instead.
  femib::types::mesh<float, 2> mesh(
      {{0, 0}, {2, 0}, {2, 1}, {1, 1}, {1, 2}, {0, 2}},
      {{0, 1, 2}, {0, 2, 3}, {0, 3, 4}, {0, 4, 5}}, {0, 1, 2, 3, 4, 5});

  femib::finite_element::finite_element<float, 2, 2> f =
      femib::finite_element::create_finite_element_P1_2d2d<float, 2, 2>();
  femib::finite_element_space::finite_element_space<float, 2, 2> s(f, mesh);

  Eigen::Matrix<float, Eigen::Dynamic, 1> xx =
      Eigen::Matrix<float, Eigen::Dynamic, 1>::Zero(s.nodes.P.size(), 1);

  CHECK_NOTHROW(s.plot(xx, 0.01));
}

TEST_CASE("interpolate is exact for a linear field at points rounding-close "
          "to mesh edges") {
  using vec = femib::types::dvec<double, 2>;
  // two triangles sharing the diagonal x + y = 1; points just across it, in
  // the second one, are close enough for find_points to give the first
  femib::types::mesh<double, 2> mesh(
      {vec(0, 0), vec(1, 0), vec(0, 1), vec(1, 1)},
      {femib::types::ditrian<2>(0, 1, 2), femib::types::ditrian<2>(1, 3, 2)});
  femib::finite_element::finite_element<double, 2, 1> f =
      femib::finite_element::create_finite_element_P1_2d1d<double, 2, 1>();
  femib::finite_element_space::finite_element_space<double, 2, 1> s(f, mesh);

  auto linear = [](const vec &p) { return 1 + p(0) + 2 * p(1); };
  Eigen::Matrix<double, Eigen::Dynamic, 1> u(s.size());
  for (int i = 0; i < s.size(); ++i) {
    u(i) = linear(s.nodes.P[i]);
  }

  std::vector<vec> points;
  for (double t : {0.1, 0.3, 0.5, 0.7, 0.9}) {
    points.push_back(vec(t, 1 - t) + 1e-12 * vec(1, 1));
  }
  std::vector<int> located = femib::mesh::find_points<double, 2>(mesh, points);
  std::vector<femib::types::dvec<double, 1>> values = s.interpolate(u, points);

  for (size_t k = 0; k < points.size(); ++k) {
    CAPTURE(k);
    CHECK_GE(located[k], 0);
    CHECK(values[k](0) == doctest::Approx(linear(points[k])).epsilon(1e-9));
  }
}
