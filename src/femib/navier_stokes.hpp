#ifndef FEMIB_NAVIER_STOKES_HPP_INCLUDED_
#define FEMIB_NAVIER_STOKES_HPP_INCLUDED_

#include "../affine/affine.hpp"
#include "../femib/femib.hpp"
#include "../finite_element_space/finite_element_space.hpp"
#include "../gauss/gauss.hpp"
#include "../types/differential_operation.hpp"
#include "stokes.hpp"
#include "stokes_t.hpp"
#include <Eigen/Dense>
#include <algorithm>
#include <optional>

namespace femib::navier_stokes {

// convection term, using Picard's iterative method
template <typename T, int d>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> assemble_convection(
    femib::finite_element_space::finite_element_space<T, d, d> &V,
    const femib::gauss::rule<T, d> &rule,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>
        &w_dofs, // TODO uh? u_at_last_step
    T reynolds   // TODO reynolds? nu!
) {
  std::vector<Eigen::Triplet<T>> BB;
  for (int n = 0; n < V.mesh.T.size(); ++n) {
    femib::types::dtrian<T, d> t = V.mesh[n];
    femib::types::dmat<T, d> Binv = femib::affine::affineBinv(t);
    femib::types::dvec<T, d> bb = femib::affine::affineb(t);
    femib::types::F<T, d, d> w;
    w.x = [&V, n, Binv, bb, &w_dofs](
              const femib::types::dvec<T, d> &x) -> femib::types::dvec<T, d> {
      femib::types::dvec<T, d> res = femib::types::dvec<T, d>::Zero();
      femib::types::dvec<T, d> x_ref = Binv * (x - bb);
      for (int k = 0; k < V.finite_element.base_functions.size(); ++k) {
        int global_k = V.nodes.get_index(k, n);
        res += w_dofs(global_k) * V.finite_element.base_functions[k].x(x_ref);
      }
      return res;
    };
    w.dx = [](const femib::types::dvec<T, d> &) -> femib::types::rmat<T, d, d> {
      return femib::types::rmat<T, d, d>::Zero();
    };
    for (int i = 0; i < V.finite_element.base_functions.size(); ++i) {
      femib::types::F<T, d, d> a =
          femib::util::base_function2real_function<T, d, d>(V, i, Binv, bb);
      for (int j = 0; j < V.finite_element.base_functions.size(); ++j) {
        femib::types::F<T, d, d> b =
            femib::util::base_function2real_function<T, d, d>(V, j, Binv, bb);
        auto fff = [&](const femib::types::dvec<T, d> &x) {
          return a.x(x).dot(convection<T, d>(w, b)(x));
        };
        T m = femib::mesh::integrate<T, d>(rule, fff, t);
        BB.push_back(Eigen::Triplet<T>(V.nodes.get_index(i, n),
                                       V.nodes.get_index(j, n), m));
      }
    }
  }
  return (T(1) / reynolds) *
         femib::util::triplets2dense(BB, V.nodes.P.size(), V.nodes.P.size());
}

// Picard iteration:
// solve Stokes once, then repeatedly add the convection contribution to the
// frozen Stokes matrix, rebuild the solvable system, and re-solve, until the
// relative change is below tol or max_picard_iters is reached.
template <typename T, int d>
Eigen::Matrix<T, Eigen::Dynamic, 1>
solve_steady(femib::stokes::stokes<T, d> &s,
             const femib::gauss::rule<T, d> &rule, T reynolds,
             int max_picard_iters, T tol) {
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A_stokes = s.A;
  Eigen::Matrix<T, Eigen::Dynamic, 1> xx = femib::stokes::solve<T, d, 1>(s);
  for (int iter = 0; iter < max_picard_iters; ++iter) {
    Eigen::Matrix<T, Eigen::Dynamic, 1> w_dofs = xx.topRows(s.V.nodes.P.size());
    s.A = A_stokes + assemble_convection<T, d>(s.V, rule, w_dofs, reynolds);
    femib::stokes::rebuild_system<T, d>(s, rule);
    Eigen::Matrix<T, Eigen::Dynamic, 1> xx_new =
        femib::stokes::solve<T, d, 1>(s);
    T rel_change =
        (xx_new - xx).norm() / std::max(xx_new.norm(), static_cast<T>(1e-8));
    xx = xx_new;
    if (rel_change < tol)
      break;
  }
  return xx;
}

// backward-Euler timestep of the time-dependent Navier-Stokes system,
// with an inner Picard loop
template <typename T, int d>
void advance(femib::stokes_t::stokes<T, d> &s,
             const femib::gauss::rule<T, d> &rule, T reynolds,
             int max_picard_iters, T tol,
             std::optional<Eigen::Matrix<T, Eigen::Dynamic, 1>>
                 extra_velocity_rhs = std::nullopt) {
  s.time += s.deltat;

  Eigen::Matrix<T, Eigen::Dynamic, 1> u_1;
  if (s.solution.size() == 0)
    u_1 = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(s.V.nodes.P.size(), 1);
  else
    u_1 = s.solution[s.solution.size() - 1].topRows(s.V.nodes.P.size());
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> DD =
      (1 / s.deltat) *
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(
          s.V.nodes.P.size(), s.V.nodes.P.size());
  Eigen::Matrix<T, Eigen::Dynamic, 1> dd = (1 / s.deltat) * u_1;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A_base = s.A + DD;

  // Picard loop
  Eigen::Matrix<T, Eigen::Dynamic, 1> xx = u_1;
  Eigen::Matrix<T, Eigen::Dynamic, 1> xx_new_full;
  for (int iter = 0; iter < max_picard_iters; ++iter) {
    s.AA.block(0, 0, s.V.nodes.P.size(), s.V.nodes.P.size()) =
        A_base + assemble_convection<T, d>(s.V, rule, xx, reynolds);

    femib::stokes_t::rebuild_system<T, d>(s, dd, extra_velocity_rhs);
    xx_new_full = femib::stokes_t::solve<T, d, 1>(s);

    Eigen::Matrix<T, Eigen::Dynamic, 1> xx_new_velocity =
        xx_new_full.topRows(s.V.nodes.P.size());
    T rel_change = (xx_new_velocity - xx).norm() /
                   std::max(xx_new_velocity.norm(), static_cast<T>(1e-8));
    xx = xx_new_velocity;
    if (iter > 0 && rel_change < tol) {
      break;
    }
  }
  s.plotV.emplace_back(s.V.plot(xx_new_full.topRows(s.V.nodes.P.size()), 0.01));
  s.plotQ.emplace_back(s.Q.plot(xx_new_full.tail(s.Q.nodes.P.size()), 0.01));
  s.solution.emplace_back(xx_new_full);
}

} // namespace femib::navier_stokes
#endif
