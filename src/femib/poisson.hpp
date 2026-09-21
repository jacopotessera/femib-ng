#ifndef FEMIB_POISSON_HPP_INCLUDED_
#define FEMIB_POISSON_HPP_INCLUDED_

#include "../affine/affine.hpp"
#include "../finite_element_space/finite_element_space.hpp"
#include "../gauss/gauss.hpp"
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <algorithm>
#include <functional>
#include <iostream>
#include <memory>
#include <spdlog/spdlog.h>
#include <vector>

#include "femib.hpp"

namespace femib::poisson {

// Solves A x = b for a sparse A. A concrete implementation owns its own
// strategy/parameters (e.g. which iterative method, tolerance); callers only
// ever see this interface.
  // TODO move to new file: sparse_solvers.hpp
template <typename T> struct sparse_solver {
  virtual ~sparse_solver() = default;
  virtual Eigen::Matrix<T, Eigen::Dynamic, 1>
  solve(const Eigen::SparseMatrix<T> &A,
        const Eigen::Matrix<T, Eigen::Dynamic, 1> &b) = 0;
};

  // TODO move to new file: sparse_solvers.hpp
template <typename T> struct bicgstab_ilut_solver : sparse_solver<T> {
  T tolerance = 1e-10;
  int max_iterations = 500;

  Eigen::Matrix<T, Eigen::Dynamic, 1>
  solve(const Eigen::SparseMatrix<T> &A,
        const Eigen::Matrix<T, Eigen::Dynamic, 1> &b) override {
    Eigen::BiCGSTAB<Eigen::SparseMatrix<T>, Eigen::IncompleteLUT<T>> solver;
    solver.setTolerance(tolerance);
    solver.setMaxIterations(max_iterations);
    solver.compute(A);

    Eigen::Matrix<T, Eigen::Dynamic, 1> x = solver.solve(b);
    if (solver.info() != Eigen::Success) {
      // TODO generic message
      SPDLOG_ERROR("femib::poisson::bicgstab_ilut_solver: BiCGSTAB+ILUT "
                   "failed to converge (Eigen::ComputationInfo = {}, "
                   "iterations = {}, estimated error = {})",
                   static_cast<int>(solver.info()), solver.iterations(),
                   solver.error());
    }
    return x;
  }
};

  // TODO move to new file: sparse_solvers.hpp
template <typename T>
std::unique_ptr<sparse_solver<T>> default_solver_factory() {
  return std::make_unique<bicgstab_ilut_solver<T>>();
}

template <typename T, int d, int e> struct poisson {

  femib::finite_element_space::finite_element_space<T, d, e> V;
  // TODO we need default here and also in the constructor?
  std::function<std::function<T(femib::types::dvec<T, d>)>(
      femib::types::F<T, d, e>)>
      force = [](femib::types::F<T, d, e>) {
        return [](femib::types::dvec<T, d>) { return T(0); };
      };

  // TODO we need default here and also in the constructor?
  std::unique_ptr<sparse_solver<T>> solver = default_solver_factory<T>();

  Eigen::Matrix<T, Eigen::Dynamic, 1> dB;
  Eigen::SparseMatrix<T> dM;
  Eigen::Matrix<T, Eigen::Dynamic, 1> dF;

  poisson() = default;
  poisson(
      femib::finite_element_space::finite_element_space<T, d, e> v,
      const femib::gauss::rule<T, d> &rule,
      std::unique_ptr<sparse_solver<T>> solver = default_solver_factory<T>())
      : V(std::move(v)), solver(std::move(solver)) { // TODO can we replace move with reference?
    std::function<T(femib::types::dvec<T, d>)> zero_boundary =
        [](const femib::types::dvec<T, d> &x) { return 0; }; // TODO we need to take this in the constructor, default 0 is ok
    std::vector<Eigen::Triplet<T>> B =
        femib::util::build_edges<T, d, e>(V, zero_boundary);

    dB = femib::util::triplets2dense<T>(B, V.nodes.P.size(), 1);
    dM = femib::util::build_stiffness_matrix(V, rule);
    dF = femib::util::build_load_vector(V, rule, force);
  }
};

  // TODO can this go to femib.hpp?
template <typename T>
femib::util::sparse_solvable_equations<T>
remove_edges(const Eigen::SparseMatrix<T> &dM,
             const Eigen::Matrix<T, Eigen::Dynamic, 1> &dF,
             const Eigen::Matrix<T, Eigen::Dynamic, 1> &dB,
             const std::vector<int> &not_edges) {

  Eigen::Matrix<T, Eigen::Dynamic, 1> ss = dM * dB;

  int n = static_cast<int>(not_edges.size());
  Eigen::Matrix<T, Eigen::Dynamic, 1> bbb(n);
  for (int i = 0; i < n; ++i) {
    bbb(i) = dF(not_edges[i]) - ss(not_edges[i]);
  }

  Eigen::SparseMatrix<T> P =
      femib::util::selection_matrix<T>(not_edges, static_cast<int>(dM.rows()));
  Eigen::SparseMatrix<T> AAA = P * dM * P.transpose();

  return {AAA, bbb};
}

  // TODO can this go to femib.hpp?
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> add_edges(

    Eigen::Matrix<T, Eigen::Dynamic, 1> xxx,
    Eigen::Matrix<T, Eigen::Dynamic, 1> dB, int rows,
    std::vector<int> not_edges, std::vector<int> nodesE) {

  Eigen::Matrix<T, Eigen::Dynamic, 1> xx;
  xx.resize(rows, 1);

  for (int i = 0; i < rows; i++) {
    xx(i, 0) = 0.0;
    auto k = std::find(not_edges.begin(), not_edges.end(), i);
    if (k != not_edges.end()) {
      xx(i, 0) = xxx(k - not_edges.begin(), 0);
    }
    auto kk = std::find(nodesE.begin(), nodesE.end(), i);
    if (kk != nodesE.end()) {
      xx(i, 0) = dB(i, 0);
    }
  }

  return xx;
}

template <typename T, int d, int e>
Eigen::Matrix<T, Eigen::Dynamic, 1> solve(const poisson<T, d, e> &poisson) {

  std::vector<int> not_edges = femib::util::build_not_edges<T, d, e>(poisson.V);

  femib::util::sparse_solvable_equations<T> solvable_equations =
      remove_edges<T>(poisson.dM, poisson.dF, poisson.dB, not_edges);

  Eigen::Matrix<T, Eigen::Dynamic, 1> x =
      poisson.solver->solve(solvable_equations.A, solvable_equations.b);

  return add_edges<T>(x, poisson.dB, poisson.V.nodes.P.size(), not_edges,
                      poisson.V.nodes.E);
}

} // namespace femib::poisson
#endif
