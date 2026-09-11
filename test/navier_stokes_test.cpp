#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../src/affine/affine.hpp"
#include "../src/femib/navier_stokes.hpp"
#include "../src/femib/stokes_t.hpp"
#include "../src/finite_element/P0_2d1d.hpp"
#include "../src/finite_element/P1+B_2d2d.hpp"
#include "../src/finite_element/P1_2d1d.hpp"
#include "../src/finite_element/P1_2d2d.hpp"
#include "../src/finite_element_space/finite_element_space.hpp"
#include "../src/gauss/gauss.hpp"
#include "../src/gauss/gauss_lagrange_2_2d.hpp"
#include "../src/types/differential_operation.hpp"
#include <cmath>
#include <doctest/doctest.h>

#include "utils.hpp"

TEST_CASE("convection((w), (u)) matches (w . grad) u for a known polynomial "
          "field, evaluated by hand") {
  femib::types::F<float, 2, 2> w;
  w.x = [](const femib::types::dvec<float, 2> &) {
    return femib::types::dvec<float, 2>(2.0f, 3.0f);
  };
  w.dx = [](const femib::types::dvec<float, 2> &)
      -> femib::types::rmat<float, 2, 2> {
    return femib::types::rmat<float, 2, 2>::Zero();
  };

  femib::types::F<float, 2, 2> u;
  u.x = [](const femib::types::dvec<float, 2> &x) {
    return femib::types::dvec<float, 2>(x(0) * x(0), x(0) * x(1));
  };
  u.dx = [](const femib::types::dvec<float, 2> &x)
      -> femib::types::rmat<float, 2, 2> {
    femib::types::rmat<float, 2, 2> J;
    J << 2 * x(0), x(1), 0, x(0);
    return J;
  };

  auto conv = convection<float, 2>(w, u);
  femib::types::dvec<float, 2> at(1.5f, -0.5f);
  femib::types::dvec<float, 2> result = conv(at);

  CHECK_LT(std::abs(result(0) - 4.0f * at(0)), 1e-4f);
  CHECK_LT(std::abs(result(1) - (2.0f * at(1) + 3.0f * at(0))), 1e-4f);
}

TEST_CASE("solve steady Navier-Stokes") {
  constexpr float PI = std::numbers::pi_v<float>;
  float eps = 0.05f;

  auto u1_exact = [eps](const float x, const float y) {
    return eps * 2.0f * PI * std::sin(PI * x) * std::sin(PI * x) *
           std::sin(PI * y) * std::cos(PI * y);
  };
  auto u2_exact = [eps](const float x, const float y) {
    return eps * -2.0f * PI * std::sin(PI * x) * std::sin(PI * y) *
           std::sin(PI * y) * std::cos(PI * x);
  };

  auto f1 = [eps](const float x, const float y) {
    return PI * eps *
           (4.0f * PI * PI * eps * std::pow(std::sin(PI * x), 3) *
                std::sin(PI * y) * std::sin(PI * y) * std::cos(PI * x) -
            (-16.0f * PI * PI * std::sin(PI * x) * std::sin(PI * x) *
                 std::sin(PI * y) +
             4.0f * PI * PI * std::sin(PI * y) + std::cos(PI * x)) *
                std::cos(PI * y));
  };
  auto f2 = [eps](const float x, const float y) {
    return PI * eps *
           (4.0f * PI * PI * eps * std::sin(PI * x) *
                std::pow(std::sin(PI * y), 3) * std::cos(PI * y) -
            16.0f * PI * PI * std::sin(PI * y) * std::sin(PI * y) *
                std::cos(PI * x) +
            std::sin(PI * y) + 4.0f * PI * PI * std::cos(PI * x)) *
           std::sin(PI * x);
  };

  int n = 20;
  femib::types::mesh<float, 2> mesh = make_unit_square_mesh(n);
  mesh.init();

  femib::gauss::rule<float, 2> rule =
      femib::gauss::create_gauss_2_2d<float, 2>();

  femib::finite_element::finite_element<float, 2, 2> f_p1_2d2d =
      femib::finite_element::create_finite_element_P1_B_2d2d<float, 2, 2>();
  femib::finite_element_space::finite_element_space v = {
      .finite_element = f_p1_2d2d, .mesh = mesh};
  v.nodes = f_p1_2d2d.build_nodes(mesh);

  femib::finite_element::finite_element<float, 2, 1> f_p1_2d1d =
      femib::finite_element::create_finite_element_P1_2d1d<float, 2, 1>();
  femib::finite_element_space::finite_element_space q = {
      .finite_element = f_p1_2d1d, .mesh = mesh};
  q.nodes = f_p1_2d1d.build_nodes(mesh);

  femib::stokes::stokes<float, 2> s;
  s.V = v;
  s.Q = q;

  // external force
  // TODO stokes interface need improvements?
  auto ggg = [f1, f2](const femib::types::F<float, 2, 2> &a) {
    return [a, f1, f2](const femib::types::dvec<float, 2> &x) {
      return f1(x(0), x(1)) * a.x(x)(0) + f2(x(0), x(1)) * a.x(x)(1);
    };
  };

  femib::util::build_diagonal_result<float> result =
      femib::util::build_diagonal<float, 2, 2>(
          s.V, rule, femib::stokes::stokes_a<float, 2>, ggg);
  s.A = femib::util::triplets2dense(result.M, s.V.nodes.P.size(),
                                    s.V.nodes.P.size());
  s.B = femib::util::triplets2dense(
      femib::util::build_non_diagonal<float, 2>(
          s.V, s.Q, rule, femib::stokes::stokes_b<float, 2>),
      s.V.nodes.P.size(), s.Q.nodes.P.size());

  s.ff = Eigen::Matrix<float, Eigen::Dynamic, 1>::Zero(
      s.V.nodes.P.size() + s.Q.nodes.P.size(), 1);
  s.ff.block(0, 0, s.V.nodes.P.size(), 1) =
      femib::util::triplets2dense(result.F, s.V.nodes.P.size(), 1);

  femib::stokes::rebuild_system<float, 2>(s, rule);

  Eigen::Matrix<float, Eigen::Dynamic, 1> xx =
      femib::navier_stokes::solve_steady<float, 2>(s, rule, /*reynolds=*/1.0f,
                                                   /*max_picard_iters=*/30,
                                                   /*tol=*/1e-5f);

  size_t size_P = mesh.P.size();
  size_t size_T = mesh.T.size();
  float max_velocity_error = 0.0f;
  for (int i = 0; i < size_P; ++i) {
    float x = mesh.P[i](0), y = mesh.P[i](1);
    max_velocity_error =
        std::max(max_velocity_error, std::abs(xx(i) - u1_exact(x, y)));
    max_velocity_error =
        std::max(max_velocity_error, std::abs(xx(size_P + i) - u2_exact(x, y)));
  }
  for (int n_ = 0; n_ < size_T; ++n_) {
    femib::types::dvec<float, 2> c =
        femib::finite_element::find_center_of<float, 2>(mesh[n_]);
    max_velocity_error =
        std::max(max_velocity_error,
                 std::abs(xx(2 * size_P + n_) - u1_exact(c(0), c(1))));
    max_velocity_error =
        std::max(max_velocity_error,
                 std::abs(xx(2 * size_P + size_T + n_) - u2_exact(c(0), c(1))));
  }
  CHECK_LT(max_velocity_error, 0.15f * eps);
}

// TODO uh?
TEST_CASE("assemble_convection matches by-hand-derived closed-form integrals "
          "on a single reference triangle") {
  femib::types::mesh<float, 2> mesh;
  mesh.P = {femib::types::dvec<float, 2>(0.0f, 0.0f),
            femib::types::dvec<float, 2>(1.0f, 0.0f),
            femib::types::dvec<float, 2>(0.0f, 1.0f)};
  mesh.T = {femib::types::ditrian<2>(0, 1, 2)};
  mesh.E = {0, 1, 2};
  mesh.init();

  femib::gauss::rule<float, 2> rule =
      femib::gauss::create_gauss_2_2d<float, 2>();

  femib::finite_element::finite_element<float, 2, 2> f_p1_2d2d =
      femib::finite_element::create_finite_element_P1_2d2d<float, 2, 2>();
  femib::finite_element_space::finite_element_space<float, 2, 2> v = {
      .finite_element = f_p1_2d2d, .mesh = mesh};
  v.nodes = f_p1_2d2d.build_nodes(mesh);

  Eigen::Matrix<float, Eigen::Dynamic, 1> w_dofs(6);
  w_dofs << 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f;

  Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> M =
      femib::navier_stokes::assemble_convection<float, 2>(v, rule, w_dofs,
                                                          /*reynolds=*/1.0f);

  REQUIRE_EQ(M.rows(), 6);
  REQUIRE_EQ(M.cols(), 6);

  float tol = 1e-4f;
  // u1 block (rows/cols 0,1,2): M(i,j) = integral(basis_i . (w.grad)basis_j)
  CHECK_LT(std::abs(M(0, 0) - (-13.0f / 12.0f)), tol);
  CHECK_LT(std::abs(M(1, 0) - (-7.0f / 6.0f)), tol);
  CHECK_LT(std::abs(M(2, 0) - (-5.0f / 4.0f)), tol);
  CHECK_LT(std::abs(M(0, 1) - (7.0f / 24.0f)), tol);
  CHECK_LT(std::abs(M(1, 1) - (1.0f / 3.0f)), tol);
  CHECK_LT(std::abs(M(2, 1) - (3.0f / 8.0f)), tol);
  CHECK_LT(std::abs(M(0, 2) - (19.0f / 24.0f)), tol);
  CHECK_LT(std::abs(M(1, 2) - (5.0f / 6.0f)), tol);
  CHECK_LT(std::abs(M(2, 2) - (7.0f / 8.0f)), tol);
  // u2 block (rows/cols 3,4,5) is numerically identical to the u1 block:
  // basis 3,4,5's shapes (1-x-y),(x),(y) in their (second) nonzero
  // component exactly mirror basis 0,1,2's, and w's two components enter
  // symmetrically in conv_component2(j=3)=-(w1+w2), conv_component2(j=4)=w1,
  // conv_component2(j=5)=w2 -- identical in form to the u1 block's j=0,1,2.
  CHECK_LT(std::abs(M(3, 3) - (-13.0f / 12.0f)), tol);
  CHECK_LT(std::abs(M(4, 3) - (-7.0f / 6.0f)), tol);
  CHECK_LT(std::abs(M(5, 3) - (-5.0f / 4.0f)), tol);
  CHECK_LT(std::abs(M(3, 4) - (7.0f / 24.0f)), tol);
  CHECK_LT(std::abs(M(4, 4) - (1.0f / 3.0f)), tol);
  CHECK_LT(std::abs(M(5, 4) - (3.0f / 8.0f)), tol);
  CHECK_LT(std::abs(M(3, 5) - (19.0f / 24.0f)), tol);
  CHECK_LT(std::abs(M(4, 5) - (5.0f / 6.0f)), tol);
  CHECK_LT(std::abs(M(5, 5) - (7.0f / 8.0f)), tol);
  // Cross blocks (u1 test against u2 trial and vice versa) are exactly
  // zero: basis functions 0,1,2 have zero second component and 3,4,5 have
  // zero first component, so their dot products with any (w.grad)basis_j
  // vanish identically regardless of w.
  for (int i = 0; i < 3; ++i) {
    for (int j = 3; j < 6; ++j) {
      CHECK_LT(std::abs(M(i, j)), tol);
      CHECK_LT(std::abs(M(j, i)), tol);
    }
  }
  // 1/reynolds scaling.
  Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> M_half_reynolds =
      femib::navier_stokes::assemble_convection<float, 2>(v, rule, w_dofs,
                                                          /*reynolds=*/2.0f);
  CHECK_LT((M_half_reynolds - 0.5f * M).norm(), tol);
}

// TODO fix this when the solver is fast, using sparse solver
TEST_CASE("testing time-dependent Navier-Stokes" * doctest::skip(true)) {
  constexpr float PI = std::numbers::pi_v<float>;
  float eps = 0.05f;

  auto u1_exact = [eps](const float x, const float y) {
    return eps * 2.0f * PI * std::sin(PI * x) * std::sin(PI * x) *
           std::sin(PI * y) * std::cos(PI * y);
  };
  auto u2_exact = [eps](const float x, const float y) {
    return eps * -2.0f * PI * std::sin(PI * x) * std::sin(PI * y) *
           std::sin(PI * y) * std::cos(PI * x);
  };
  auto f1 = [eps](const float x, const float y) {
    return PI * eps *
           (4.0f * PI * PI * eps * std::pow(std::sin(PI * x), 3) *
                std::sin(PI * y) * std::sin(PI * y) * std::cos(PI * x) -
            (-16.0f * PI * PI * std::sin(PI * x) * std::sin(PI * x) *
                 std::sin(PI * y) +
             4.0f * PI * PI * std::sin(PI * y) + std::cos(PI * x)) *
                std::cos(PI * y));
  };
  auto f2 = [eps](const float x, const float y) {
    return PI * eps *
           (4.0f * PI * PI * eps * std::sin(PI * x) *
                std::pow(std::sin(PI * y), 3) * std::cos(PI * y) -
            16.0f * PI * PI * std::sin(PI * y) * std::sin(PI * y) *
                std::cos(PI * x) +
            std::sin(PI * y) + 4.0f * PI * PI * std::cos(PI * x)) *
           std::sin(PI * x);
  };

  int n = 20;
  femib::types::mesh<float, 2> mesh = make_unit_square_mesh(n);
  mesh.init();
  femib::gauss::rule<float, 2> rule =
      femib::gauss::create_gauss_2_2d<float, 2>();

  femib::finite_element::finite_element<float, 2, 2> f_p1_2d2d =
      femib::finite_element::create_finite_element_P1_B_2d2d<float, 2, 2>();
  femib::finite_element_space::finite_element_space<float, 2, 2> v = {
      .finite_element = f_p1_2d2d, .mesh = mesh};
  v.nodes = f_p1_2d2d.build_nodes(mesh);

  femib::finite_element::finite_element<float, 2, 1> f_p0_2d1d =
      femib::finite_element::create_finite_element_P0_2d1d<float, 2, 1>();
  femib::finite_element_space::finite_element_space<float, 2, 1> q = {
      .finite_element = f_p0_2d1d, .mesh = mesh};
  q.nodes = f_p0_2d1d.build_nodes(mesh);

  femib::stokes_t::stokes<float, 2> s;
  s.V = v;
  s.Q = q;
  s.deltat = 0.02f;
  s.force = [f1, f2](const femib::types::dvec<float, 2> &x,
                     float) -> femib::types::dvec<float, 2> {
    return femib::types::dvec<float, 2>(f1(x(0), x(1)), f2(x(0), x(1)));
  };
  femib::stokes_t::init<float, 2>(s, rule);

  size_t size_P = mesh.P.size();
  size_t size_T = mesh.T.size();
  auto velocity_error = [&](const Eigen::Matrix<float, Eigen::Dynamic, 1> &xx) {
    float max_error = 0.0f;
    for (int i = 0; i < size_P; ++i) {
      float x = mesh.P[i](0), y = mesh.P[i](1);
      max_error = std::max(max_error, std::abs(xx(i) - u1_exact(x, y)));
      max_error =
          std::max(max_error, std::abs(xx(size_P + i) - u2_exact(x, y)));
    }
    for (int n_ = 0; n_ < size_T; ++n_) {
      femib::types::dvec<float, 2> c =
          femib::finite_element::find_center_of<float, 2>(mesh[n_]);
      max_error = std::max(
          max_error, std::abs(xx(2 * size_P + n_) - u1_exact(c(0), c(1))));
      max_error = std::max(max_error, std::abs(xx(2 * size_P + size_T + n_) -
                                               u2_exact(c(0), c(1))));
    }
    return max_error;
  };

  int n_steps = 30;
  float max_velocity_error_first = -1.0f;
  float max_velocity_error_last = -1.0f;
  for (int step = 0; step < n_steps; ++step) {
    femib::navier_stokes::advance<float, 2>(s, rule, 1.0f, 3, 1e-5f);
    float e = velocity_error(s.solution.back());
    if (step == 0)
      max_velocity_error_first = e;
    max_velocity_error_last = e;
  }

  CHECK_LT(max_velocity_error_last, max_velocity_error_first);
  CHECK_LT(max_velocity_error_last, 0.155f);
}

// TODO uh?
TEST_CASE("femib::navier_stokes::advance wires assemble_convection's "
          "contribution into the Picard-updated velocity block of s.AA") {
  // The relaxation test above cannot discriminate correct convection wiring
  // from convection being silently dropped or subtly miswired: for that
  // eps=0.05 manufactured flow, assemble_convection's own contribution is
  // ~5 orders of magnitude smaller than the linear Stokes operator, so that
  // test's outcome is essentially insensitive to whether advance()'s
  // Picard loop actually adds assemble_convection's output into s.AA
  // before each solve. This test targets exactly that composition/wiring
  // question directly -- distinct from assemble_convection's own internal
  // row/column/indexing correctness, which the closed-form
  // reference-triangle test above already proves in isolation.
  //
  // Small mesh (n=6) and a large constant forcing -- this test only needs
  // a genuinely nonzero advecting field to check the wiring against, not
  // numerical accuracy.
  int n = 6;
  femib::types::mesh<float, 2> mesh = make_unit_square_mesh(n);
  mesh.init();
  femib::gauss::rule<float, 2> rule = femib::gauss::create_gauss_2_2d<float, 2>();

  femib::finite_element::finite_element<float, 2, 2> f_p1_2d2d =
      femib::finite_element::create_finite_element_P1_B_2d2d<float, 2, 2>();
  femib::finite_element_space::finite_element_space<float, 2, 2> v = {
      .finite_element = f_p1_2d2d, .mesh = mesh};
  v.nodes = f_p1_2d2d.build_nodes(mesh);

  femib::finite_element::finite_element<float, 2, 1> f_p0_2d1d =
      femib::finite_element::create_finite_element_P0_2d1d<float, 2, 1>();
  femib::finite_element_space::finite_element_space<float, 2, 1> q = {
      .finite_element = f_p0_2d1d, .mesh = mesh};
  q.nodes = f_p0_2d1d.build_nodes(mesh);

  femib::stokes_t::stokes<float, 2> s;
  s.V = v;
  s.Q = q;
  s.deltat = 0.02f;
  // A large constant forcing (rather than the struct's default (1,1)) --
  // needed because assemble_convection's contribution scales LINEARLY with
  // the advecting field's magnitude, while A_base (dominated by the
  // backward-Euler mass term (1/deltat)*I) does NOT depend on velocity at
  // all -- so no amount of scaling makes convection comparable to A_base's
  // own scale unless the advecting field itself is large. 5000x the
  // default forcing reliably clears this test's own wiring-check tolerance
  // with comfortable margin.
  s.force = [](const femib::types::dvec<float, 2> &, float)
      -> femib::types::dvec<float, 2> {
    return femib::types::dvec<float, 2>(5000.0f, 5000.0f);
  };
  femib::stokes_t::init<float, 2>(s, rule);

  size_t rowsV = s.V.nodes.P.size();
  float reynolds = 1.0f;

  // Step 1: establishes a nonzero previous-timestep velocity, starting from
  // rest -- its own convection input (u_1=0) never matters here, only its
  // OUTPUT (used as step 2's advecting field) does. max_picard_iters=1
  // keeps this deterministic: the Picard loop's single iteration always
  // uses xx=u_1, so the field advance() actually convected against is
  // exactly known to this test, not merely assumed.
  femib::navier_stokes::advance<float, 2>(s, rule, reynolds,
                                          /*max_picard_iters=*/1,
                                          /*tol=*/1e-5f);
  REQUIRE_EQ(s.solution.size(), 1);
  Eigen::Matrix<float, Eigen::Dynamic, 1> w_dofs_expected =
      s.solution.back().topRows(rowsV);
  // Sanity: step 1's output velocity must be genuinely nonzero, otherwise
  // step 2's convection contribution would be trivially zero regardless of
  // whether the wiring is correct, making this test vacuous.
  REQUIRE_GT(w_dofs_expected.norm(), 1e-6f);

  // s.A is never modified by advance() (only s.AA's block is), and
  // s.deltat is fixed, so this is exactly the A_base advance() itself
  // computes internally for step 2 below.
  Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> A_base_expected =
      s.A + (1.0f / s.deltat) *
                Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>::Identity(
                    rowsV, rowsV);

  // Step 2 is the call under test. Its Picard loop's single iteration
  // (max_picard_iters=1) must use w_dofs_expected (step 1's output,
  // captured above) as the advecting field.
  femib::navier_stokes::advance<float, 2>(s, rule, reynolds,
                                          /*max_picard_iters=*/1,
                                          /*tol=*/1e-5f);

  Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> conv_expected =
      femib::navier_stokes::assemble_convection<float, 2>(
          s.V, rule, w_dofs_expected, reynolds);
  // Sanity: the expected convection contribution must clear this test's
  // own wiring-check tolerance (below) with comfortable margin -- otherwise
  // this check could not tell correct wiring apart from broken/missing
  // wiring either.
  REQUIRE_GT(conv_expected.norm(), 0.02f);

  // s.AA is not touched again after advance()'s Picard loop finishes, so
  // its velocity block still holds exactly what the loop's one iteration
  // computed: this is the actual assertion under test.
  Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> AA_velocity_block =
      s.AA.block(0, 0, rowsV, rowsV);
  CHECK_LT((AA_velocity_block - (A_base_expected + conv_expected)).norm(),
           1e-2f);
}
