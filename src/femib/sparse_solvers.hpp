#ifndef FEMIB_SPARSE_SOLVERS_HPP_INCLUDED_
#define FEMIB_SPARSE_SOLVERS_HPP_INCLUDED_

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <memory>
#include <spdlog/spdlog.h>

namespace femib::util {

// Solves A x = b for a sparse A. A concrete implementation owns its own
// strategy/parameters (e.g. which iterative method, tolerance)
template <typename T> struct sparse_solver {
  virtual ~sparse_solver() = default;
  virtual Eigen::Matrix<T, Eigen::Dynamic, 1>
  solve(const Eigen::SparseMatrix<T> &A,
        const Eigen::Matrix<T, Eigen::Dynamic, 1> &b) = 0;
};

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
      SPDLOG_ERROR("femib::util::bicgstab_ilut_solver: BiCGSTAB+ILUT failed "
                   "to converge (Eigen::ComputationInfo = {}, "
                   "iterations = {}, estimated error = {})",
                   static_cast<int>(solver.info()), solver.iterations(),
                   solver.error());
    }
    return x;
  }
};

template <typename T>
std::unique_ptr<sparse_solver<T>> default_solver_factory() {
  return std::make_unique<bicgstab_ilut_solver<T>>();
}

} // namespace femib::util
#endif
