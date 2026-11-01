#ifndef FEMIB_POISSON_HPP_INCLUDED_
#define FEMIB_POISSON_HPP_INCLUDED_

#include "../affine/affine.hpp"
#include "../finite_element_space/finite_element_space.hpp"
#include "../gauss/gauss.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <functional>
#include <memory>
#include <vector>

#include "femib.hpp"
#include "sparse_solvers.hpp"

namespace femib::poisson {

template <typename T, int d, int e>
std::function<
    std::function<T(femib::types::dvec<T, d>)>(femib::types::F<T, d, e>)>
default_force() {
  return [](femib::types::F<T, d, e>) {
    return [](femib::types::dvec<T, d>) { return T(0); };
  };
}

template <typename T, int d>
std::function<T(femib::types::dvec<T, d>)> default_boundary() {
  return [](const femib::types::dvec<T, d> &) { return T(0); };
}

template <typename T, int d, int e> struct poisson {

  femib::finite_element_space::finite_element_space<T, d, e> solution_space;

  std::function<std::function<T(femib::types::dvec<T, d>)>(
      femib::types::F<T, d, e>)>
      body_force;

  std::unique_ptr<femib::util::sparse_solver<T>> sparse_solver;

  Eigen::Matrix<T, Eigen::Dynamic, 1> boundary_values_on_nodes;
  Eigen::SparseMatrix<T> stiffness_matrix;
  Eigen::Matrix<T, Eigen::Dynamic, 1> load_vector;

  poisson(femib::finite_element_space::finite_element_space<T, d, e>
              solution_space_,
          const femib::gauss::rule<T, d> &integration_rule,
          std::function<T(femib::types::dvec<T, d>)> boundary =
              default_boundary<T, d>(),
          std::function<std::function<T(femib::types::dvec<T, d>)>(
              femib::types::F<T, d, e>)>
              body_force_ = default_force<T, d, e>(),
          std::unique_ptr<femib::util::sparse_solver<T>> sparse_solver_ =
              femib::util::default_sparse_solver_factory<T>())
      : solution_space(std::move(solution_space_)),
        body_force(std::move(body_force_)),
        sparse_solver(std::move(sparse_solver_)) {
    std::vector<Eigen::Triplet<T>> boundary_values_triplets =
        femib::util::build_boundary_values_on_nodes<T, d, e>(solution_space,
                                                             boundary);

    boundary_values_on_nodes = femib::util::triplets2dense<T>(
        boundary_values_triplets, solution_space.nodes.P.size(),
        1); // TODO solution_space.nodes.P.size() ?
    stiffness_matrix =
        femib::util::build_stiffness_matrix(solution_space, integration_rule);
    load_vector = femib::util::build_load_vector(solution_space,
                                                 integration_rule, body_force);
  }
};

template <typename T, int d, int e>
Eigen::Matrix<T, Eigen::Dynamic, 1> solve(const poisson<T, d, e> &poisson) {

  std::vector<int> not_edge_nodes =
      femib::util::find_not_edge_nodes<T, d, e>(poisson.solution_space);

  femib::util::sparse_linear_system<T> reduced_system =
      femib::util::build_reduced_system<T>(
          poisson.stiffness_matrix, poisson.load_vector,
          poisson.boundary_values_on_nodes, not_edge_nodes);

  Eigen::Matrix<T, Eigen::Dynamic, 1> x =
      poisson.sparse_solver->solve(reduced_system.A, reduced_system.b);

  return femib::util::build_full_solution<T>(
      x, poisson.boundary_values_on_nodes,
      poisson.solution_space.nodes.P.size(), not_edge_nodes,
      poisson.solution_space.nodes.E);
}

} // namespace femib::poisson
#endif
