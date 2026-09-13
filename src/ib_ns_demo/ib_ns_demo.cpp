// Standalone demo: an elliptical Lagrangian ring relaxing toward a circle,
// immersed in a Navier-Stokes fluid, via femib::ib::advance_navier_stokes.
// Dumps the ring's point positions, the fluid velocity field, and the
// pressure field at checkpoint intervals, to CSV files for
// plotting/animation. Not part of the test suite -- a one-off
// visualization driver.

#include "../femib/stokes_t.hpp"
#include "../finite_element/P0_2d1d.hpp"
#include "../finite_element/P1+B_2d2d.hpp"
#include "../gauss/gauss_lagrange_2_2d.hpp"
#include "../ib/ib.hpp"

#include <cstdio>
#include <fstream>
#include <iostream>

namespace {
femib::types::mesh<double, 2> make_unit_square_mesh(int n) {
  femib::types::mesh<double, 2> mesh;
  int side = n + 1;
  auto idx = [side](int i, int j) { return i * side + j; };
  for (int i = 0; i < side; ++i)
    for (int j = 0; j < side; ++j) {
      double x = static_cast<double>(i) / n;
      double y = static_cast<double>(j) / n;
      mesh.P.push_back(femib::types::dvec<double, 2>(x, y));
    }
  // Checkerboard diagonal pattern (alternate by (i+j) parity), NOT the
  // fixed-diagonal-everywhere pattern this helper's own test/ib_test.cpp
  // counterpart uses. A single fixed diagonal direction for every cell
  // breaks exact mirror symmetry of the discrete fluid solver under
  // x->1-x/y->1-y reflection (a real structured-mesh artifact, not a
  // physics/forcing issue) -- for even n, alternating by (i+j)%2 makes each
  // cell's diagonal swap to its mirrored orientation under reflection,
  // restoring exact discrete symmetry for a symmetric initial condition.
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) {
      int p00 = idx(i, j), p10 = idx(i + 1, j), p01 = idx(i, j + 1),
          p11 = idx(i + 1, j + 1);
      if ((i + j) % 2 == 0) {
        mesh.T.push_back(femib::types::ditrian<2>(p00, p10, p11));
        mesh.T.push_back(femib::types::ditrian<2>(p00, p11, p01));
      } else {
        mesh.T.push_back(femib::types::ditrian<2>(p00, p10, p01));
        mesh.T.push_back(femib::types::ditrian<2>(p10, p11, p01));
      }
    }
  for (int i = 0; i < side; ++i)
    for (int j = 0; j < side; ++j)
      if (i == 0 || i == n || j == 0 || j == n)
        mesh.E.push_back(idx(i, j));
  return mesh;
}

// Shoelace formula: signed area of the closed polygon X[0..n). The ring is
// stored in a consistent (originally counter-clockwise, from build_ring's
// increasing-theta construction) orientation, so this is positive; std::abs
// guards against a sign flip if the orientation is ever inverted.
double polygon_area(const std::vector<femib::types::dvec<double, 2>> &X) {
  int n = (int)X.size();
  double A = 0.0;
  for (int k = 0; k < n; ++k) {
    int kp1 = (k + 1) % n;
    A += X[k](0) * X[kp1](1) - X[kp1](0) * X[k](1);
  }
  return std::abs(A) * 0.5;
}
} // namespace

int main(int argc, char **argv) {
  std::string ring_path = argc > 1 ? argv[1] : "/tmp/ib_ns_demo_ring.csv";
  std::string field_path = argc > 2 ? argv[2] : "/tmp/ib_ns_demo_field.csv";
  std::string pressure_path =
      argc > 3 ? argv[3] : "/tmp/ib_ns_demo_pressure.csv";
  int n_mesh = 32;           // finer Eulerian mesh -- tractable at full run
                             // length thanks to the sparse-solver fix (this
                             // used to take ~20 min/step; now ~1s/step)
  int n_ring = 96;           // keep Lagrangian point spacing well under half
                             // the mesh spacing
  double radius = 0.15;
  double k_spring = 16000.0; // pushed higher, for a visible rebound
  double viscosity = 0.3;    // multiplies s.A (see below) -- this codebase has no
                             // wired-up viscosity coefficient (mu=2 is baked into
                             // dpi(symm(u),symm(v)) with no adjustable multiplier,
                             // a known gap), so this scales the assembled viscous
                             // stiffness matrix directly, locally, in the demo.
  double deltat = 0.0002;
  double reynolds = 5.0;
  int max_picard_iters = 5;
  double tol = 1e-4;
  int n_steps = 360;         // enough to capture the full bounce + settle
  int field_every = 10;      // checkpoint interval for fluid-field snapshots

  femib::gauss::rule<double, 2> rule =
      femib::gauss::create_gauss_2_2d<double, 2>();
  femib::types::mesh<double, 2> mesh = make_unit_square_mesh(n_mesh);
  mesh.init();

  femib::finite_element::finite_element<double, 2, 2> f_p1_2d2d =
      femib::finite_element::create_finite_element_P1_B_2d2d<double, 2, 2>();
  femib::finite_element_space::finite_element_space<double, 2, 2> v = {
      .finite_element = f_p1_2d2d, .mesh = mesh};
  v.nodes = f_p1_2d2d.build_nodes(mesh);

  femib::finite_element::finite_element<double, 2, 1> f_p0_2d1d =
      femib::finite_element::create_finite_element_P0_2d1d<double, 2, 1>();
  femib::finite_element_space::finite_element_space<double, 2, 1> q = {
      .finite_element = f_p0_2d1d, .mesh = mesh};
  q.nodes = f_p0_2d1d.build_nodes(mesh);

  femib::ib::ib_problem<double, 2> p;
  p.fluid.V = v;
  p.fluid.Q = q;
  p.fluid.deltat = deltat;
  p.fluid.force = [](const femib::types::dvec<double, 2> &, double)
      -> femib::types::dvec<double, 2> {
    return femib::types::dvec<double, 2>::Zero();
  };
  p.structure = femib::ib::build_ring<double, 2>(
      femib::types::dvec<double, 2>(0.5, 0.5), radius, n_ring, k_spring);

  // This codebase's ring is already a zero-rest-length spring (force =
  // k * (X[i_next] - X[i]), always pulling inward -- never reaches a
  // zero-force state on its own except at a single point), so unlike a
  // general nonzero-rest-length model there is no separate rest_length
  // field to zero out here. What stops the ring from collapsing to a
  // point is the fluid's incompressibility: the ring is a material curve
  // advected by a divergence-free velocity field, so the AREA it encloses
  // is conserved by the flow. Among all closed curves enclosing a FIXED
  // area, the circle uniquely minimizes total perimeter (the isoperimetric
  // inequality) -- so minimizing tension energy subject to fixed enclosed
  // area has a unique global minimizer: the circle. This is a real
  // physical model (a soap-film/surface-tension membrane), not just a
  // numerical trick.

  femib::ib::init<double, 2>(p, rule);

  // Scale the assembled viscous stiffness matrix directly -- init() has
  // already built s.A (the dpi(symm(u),symm(v)) viscous bilinear form) and
  // s.M (the mass matrix), but every advance() call re-derives s.AA's
  // velocity block from s.A fresh each timestep, so scaling s.A here, once,
  // before the timestepping loop starts, is sufficient -- every subsequent
  // call picks up the scaled value. s.B (the divergence/incompressibility
  // coupling) is untouched, so the fluid stays exactly incompressible;
  // only the viscous dissipation is weakened.
  p.fluid.A *= viscosity;

  // Perturb the circle into a dramatic ellipse (stretch x by 1.8, compress y
  // by 1/1.8, preserving enclosed area to first order), placing points at
  // EQUAL ARC LENGTH along the ellipse rather than equal angle (avoids a
  // spurious local point-spacing energy that would otherwise dominate and
  // mask the real ellipse->circle relaxation). Center is EXACTLY (0.5,0.5)
  // (matching build_ring's own center above, and the domain's symmetry
  // point) so the whole problem stays symmetric under x->1-x/y->1-y --
  // together with the checkerboard mesh above, this makes the discrete
  // problem itself symmetric, not just the continuous one. k_spring/dS
  // (material properties of the rest configuration) are untouched; only
  // the perturbed positions change.
  femib::types::dvec<double, 2> center(0.5, 0.5);
  double a_axis = 1.8 * radius;
  double b_axis = radius / 1.8;
  int n_dense = 20000;
  std::vector<double> cum_arc(n_dense + 1, 0.0);
  std::vector<double> theta_dense(n_dense + 1);
  const double PI = 3.14159265358979;
  for (int i = 0; i <= n_dense; ++i) {
    theta_dense[i] = 2.0 * PI * i / n_dense;
  }
  for (int i = 1; i <= n_dense; ++i) {
    auto speed = [&](double th) {
      double dx = -a_axis * std::sin(th);
      double dy = b_axis * std::cos(th);
      return std::sqrt(dx * dx + dy * dy);
    };
    double mid = 0.5 * (theta_dense[i - 1] + theta_dense[i]);
    cum_arc[i] = cum_arc[i - 1] + speed(mid) * (theta_dense[i] - theta_dense[i - 1]);
  }
  double total_arc = cum_arc[n_dense];
  int n_ring_actual = (int)p.structure.X.size();
  int j = 0;
  for (int k = 0; k < n_ring_actual; ++k) {
    // Half-slot phase shift (k+0.5 instead of k): with an exactly-centered
    // ring, an integer-k parametrization puts points exactly at theta=0,
    // pi/2, pi, 3pi/2 -- exactly on the x=0.5/y=0.5 mesh lines, where
    // point-location can fail and silently return zero velocity, freezing
    // those points. The half-slot shift avoids every cardinal angle while
    // preserving the ellipse's own point-symmetry about its center (a
    // uniform phase shift of all points leaves the material curve itself
    // unchanged).
    double target = total_arc * (k + 0.5) / n_ring_actual;
    while (j < n_dense && cum_arc[j + 1] < target) {
      ++j;
    }
    double theta = theta_dense[j];
    p.structure.X[k](0) = center(0) + a_axis * std::cos(theta);
    p.structure.X[k](1) = center(1) + b_axis * std::sin(theta);
  }

  double area0 = polygon_area(p.structure.X);

  std::ofstream ring_out(ring_path);
  ring_out << "step,point_index,x,y\n";
  auto dump_ring = [&](int step) {
    for (size_t k = 0; k < p.structure.X.size(); ++k) {
      ring_out << step << "," << k << "," << p.structure.X[k](0) << ","
                << p.structure.X[k](1) << "\n";
    }
  };
  dump_ring(0);

  std::ofstream field_out(field_path);
  field_out << "step,x,y,u,v\n";
  auto dump_field = [&](int step) {
    // p.fluid.plotV is populated on every femib::navier_stokes::advance
    // call (see its own s.V.plot(...) call), so this just reads the most
    // recent snapshot rather than re-sampling.
    if (p.fluid.plotV.empty())
      return;
    for (const auto &[position, velocity] : p.fluid.plotV.back()) {
      field_out << step << "," << position(0) << "," << position(1) << ","
                 << velocity(0) << "," << velocity(1) << "\n";
    }
  };

  std::ofstream pressure_out(pressure_path);
  pressure_out << "step,x,y,pressure\n";
  auto dump_pressure = [&](int step) {
    // p.fluid.plotQ is populated on every femib::navier_stokes::advance
    // call (see its own s.Q.plot(...) call) alongside plotV above.
    if (p.fluid.plotQ.empty())
      return;
    for (const auto &[position, pressure] : p.fluid.plotQ.back()) {
      pressure_out << step << "," << position(0) << "," << position(1) << ","
                    << pressure(0) << "\n";
    }
  };

  auto aspect_of = [&]() {
    double xmin = 1e9, xmax = -1e9, ymin = 1e9, ymax = -1e9;
    for (const auto &x : p.structure.X) {
      xmin = std::min(xmin, x(0));
      xmax = std::max(xmax, x(0));
      ymin = std::min(ymin, x(1));
      ymax = std::max(ymax, x(1));
    }
    return (xmax - xmin) / (ymax - ymin);
  };

  for (int step = 1; step <= n_steps; ++step) {
    femib::ib::advance_navier_stokes<double, 2>(p, rule, reynolds,
                                               max_picard_iters, tol);
    dump_ring(step);
    if (step % field_every == 0 || step == 1) {
      dump_field(step);
      dump_pressure(step);
    }
    if (step % 20 == 0 || step == 1) {
      double E = femib::ib::elastic_energy<double, 2>(p.structure);
      double area = polygon_area(p.structure.X);
      std::cerr << "step " << step << "/" << n_steps << "  elastic_energy="
                << E << "  aspect=" << aspect_of() << "  area=" << area
                << "  area/area0=" << (area / area0) << std::endl;
    }
    if (step % field_every == 0) {
      // Flush periodically (not every step -- I/O overhead) so a reader
      // polling the CSVs mid-run (e.g. to regenerate a progress GIF) always
      // sees complete, up-to-date data rather than whatever is still
      // sitting in ofstream's internal buffer.
      ring_out.flush();
      field_out.flush();
      pressure_out.flush();
    }
  }

  ring_out.close();
  field_out.close();
  pressure_out.close();
  std::cerr << "Wrote " << ring_path << ", " << field_path << ", and "
            << pressure_path << std::endl;
  return 0;
}
