// Standalone demo: an elliptical Lagrangian ring relaxing toward a circle,
// immersed in a Navier-Stokes fluid, via femib::ib::advance_navier_stokes.
// Persists the ring's point positions, the fluid velocity field, and the
// pressure field at every step (structure) / checkpoint intervals (fields)
// via this codebase's normal HDF5 persistence (femib::write::save_sim /
// save_plot_data), the same one test/ib_test.cpp uses. Not part of the test
// suite -- a one-off visualization driver; see plot/plot_simulation.py to
// render it.

#include "../femib/stokes_t.hpp"
#include "../finite_element/P0_2d1d.hpp"
#include "../finite_element/P1+B_2d2d.hpp"
#include "../gauss/gauss_lagrange_2_2d.hpp"
#include "../ib/ib.hpp"
#include "../write/write.hpp"

#include <chrono>
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
  std::string path = argc > 1 ? argv[1] : "/tmp/ib_ns_demo.h5";
  // n_mesh=32 (2*n_mesh^2 = 2048 triangles): femib::cuda::parallel_accurate
  // used to launch one CUDA thread per mesh triangle within a SINGLE block
  // (blockDim.x = size_T), which silently broke for any mesh with more
  // than ~1024 triangles (every current CUDA architecture's threads-per-
  // block cap) -- past that, find_points quietly returned "not found" for
  // every query point, and both the plotted velocity field and
  // interpolate_velocity (which advects the ring) read back as zero,
  // freezing the simulation with no error. Fixed in src/cuda/cuda.cu
  // (parallel_accurate_kernel now chunks the triangle dimension across a
  // 2D grid instead of requiring it to fit in one block), so this mesh
  // density is safe again -- finer than the original n_mesh=16, tractable
  // at full run length thanks to the sparse-solver fix (this used to take
  // ~20 min/step; now a fraction of a second per step).
  int n_mesh = 32;
  int n_ring = 2048; // keep Lagrangian point spacing well under half
                     // the mesh spacing
  double radius = 0.15;
  double k_spring = 60.0;
  double rho = 1.0;
  double mu = 0.3;
  double deltat = 0.0002;
  int max_picard_iters = 5;
  double tol = 1e-4;
  int n_steps = 480;

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
  p.fluid.rho = rho;
  p.fluid.mu = mu;
  p.fluid.force = [](const femib::types::dvec<double, 2> &,
                     double) -> femib::types::dvec<double, 2> {
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
    cum_arc[i] =
        cum_arc[i - 1] + speed(mid) * (theta_dense[i] - theta_dense[i - 1]);
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

  femib::write::save_sim(path, "ib_ns_demo", "ring",
                         femib::write::mode::overwrite);
  femib::write::save_metadata<double>(path, "rho", rho);
  femib::write::save_metadata<double>(path, "mu", mu);
  femib::write::save_metadata<double>(path, "deltat", deltat);
  femib::write::save_metadata<double>(path, "k_spring", k_spring);

  // Persists the structure's position AND the velocity/pressure fields at
  // EVERY step (so the animation has a full field for every frame, not just
  // every field_every-th one) -- femib::write::write_if_present (write.cpp)
  // silently skips any empty vector, so a timestep group with no "x"/"u"/"q"
  // datasets this step (only true at step 0, before the first solve) is
  // expected, not an error; see plot/plot_utils.py's calc_plot_data, which
  // already handles that.
  //
  // Field interpolation (velocity/pressure onto a plotting grid) used to run
  // unconditionally inside advance() at delta=0.01 (~10k points) -- the
  // actual dominant cost of a timestep once the sparse solve was fixed (see
  // the timing report). Computing it here, on demand, at every step (instead
  // of gating it behind field_every) needs a much coarser grid to stay cheap
  // -- delta=0.04 (~26x26 = 676 points) is dense enough to show the flow
  // pattern while costing ~15x fewer point-location queries than 0.01 did.
  auto dump_step = [&](int step) {
    femib::write::plot_data<double, 2> data;
    data.time = step;
    if (!p.fluid.solution.empty()) {
      for (const auto &[position, velocity] : p.fluid.plot_velocity(0.04)) {
        data.x.push_back(position);
        data.u.push_back(velocity);
      }
      for (const auto &[position, pressure] : p.fluid.plot_pressure(0.04)) {
        data.q.push_back(pressure);
      }
    }
    for (const auto &val : p.structure.X) {
      data.X.push_back(val);
    }
    femib::write::save_plot_data<double, 2>(path, data);
  };
  dump_step(0);

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

  auto t_start = std::chrono::steady_clock::now();
  for (int step = 1; step <= n_steps; ++step) {
    femib::ib::advance_navier_stokes<double, 2>(p, rule, max_picard_iters, tol);
    dump_step(step);
    if (step % 20 == 0 || step == 1) {
      double E = femib::ib::elastic_energy<double, 2>(p.structure);
      double area = polygon_area(p.structure.X);
      double elapsed_s = std::chrono::duration<double>(
                             std::chrono::steady_clock::now() - t_start)
                             .count();
      std::cerr << "step " << step << "/" << n_steps << "  elastic_energy=" << E
                << "  aspect=" << aspect_of() << "  area=" << area
                << "  area/area0=" << (area / area0)
                << "  elapsed=" << elapsed_s << "s" << std::endl;
    }
  }

  std::cerr << "Wrote " << path << std::endl;
  return 0;
}
