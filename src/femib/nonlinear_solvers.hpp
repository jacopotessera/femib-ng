#ifndef FEMIB_NONLINEAR_SOLVERS_HPP_INCLUDED_
#define FEMIB_NONLINEAR_SOLVERS_HPP_INCLUDED_

#include <Eigen/Dense>
#include <algorithm>
#include <functional>
#include <memory>
#include <stdexcept>

namespace femib::util {

template <typename T> struct nonlinear_solver {
  using vector_t = Eigen::Matrix<T, Eigen::Dynamic, 1>;
  virtual ~nonlinear_solver() = default;
  virtual vector_t
  solve(const vector_t &initial_guess,
        const std::function<vector_t(const vector_t &)> &step,
        const std::function<vector_t(const vector_t &)> &next_guess) = 0;
};

// Picard iterative method: solves the linearized problem at (guess)
template <typename T> struct picard_solver : nonlinear_solver<T> {
  using typename nonlinear_solver<T>::vector_t;
  int max_iters;
  T tol;

  explicit picard_solver(int max_iters = 20, T tol = T(1e-6))
      : max_iters(max_iters), tol(tol) {}

  vector_t
  solve(const vector_t &initial_guess,
        const std::function<vector_t(const vector_t &)> &step,
        const std::function<vector_t(const vector_t &)> &next_guess) override {
    if (max_iters <= 0) {
      throw std::invalid_argument(
          "femib::util::picard_solver: max_iters must be >= 1");
    }
    vector_t guess = initial_guess;
    vector_t full;
    for (int iter = 0; iter < max_iters; ++iter) {
      full = step(guess);
      vector_t new_guess = next_guess(full);
      T rel_change = (new_guess - guess).norm() /
                     std::max(new_guess.norm(), static_cast<T>(1e-8));
      guess = new_guess;
      if (iter > 0 && rel_change < tol) {
        break;
      }
    }
    return full;
  }
};

template <typename T>
std::unique_ptr<nonlinear_solver<T>> default_nonlinear_solver_factory() {
  return std::make_unique<picard_solver<T>>();
}

} // namespace femib::util
#endif
