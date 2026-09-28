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

  femib::finite_element_space::finite_element_space<T, d, e> V;

  std::function<std::function<T(femib::types::dvec<T, d>)>(
      femib::types::F<T, d, e>)>
      force = default_force<T, d, e>();

  std::unique_ptr<femib::util::sparse_solver<T>> solver =
      femib::util::default_solver_factory<T>();

  Eigen::Matrix<T, Eigen::Dynamic, 1> edge_values;
  Eigen::SparseMatrix<T> stiffness_matrix;
  Eigen::Matrix<T, Eigen::Dynamic, 1> load_vector;

  poisson() = default;
  poisson(femib::finite_element_space::finite_element_space<T, d, e> v,
          const femib::gauss::rule<T, d> &integration_rule,
          std::function<T(femib::types::dvec<T, d>)> boundary =
              default_boundary<T, d>(),
          std::function<std::function<T(femib::types::dvec<T, d>)>(
              femib::types::F<T, d, e>)>
              force = default_force<T, d, e>(),
          std::unique_ptr<femib::util::sparse_solver<T>> solver =
              femib::util::default_solver_factory<T>())
      : V(std::move(v)), force(std::move(force)), solver(std::move(solver)) {
    std::vector<Eigen::Triplet<T>> edge_values_triplets =
        femib::util::build_edges<T, d, e>(V, boundary);

    edge_values = femib::util::triplets2dense<T>(
        edge_values_triplets, V.nodes.P.size(), 1); // TODO V.nodes.P.size() ?
    stiffness_matrix = femib::util::build_stiffness_matrix(V, integration_rule);
    load_vector =
        femib::util::build_load_vector(V, integration_rule, this->force);
  }
};

template <typename T, int d, int e>
Eigen::Matrix<T, Eigen::Dynamic, 1> solve(const poisson<T, d, e> &poisson) {

  std::vector<int> not_edges = femib::util::build_not_edges<T, d, e>(poisson.V);

  femib::util::sparse_solvable_equations<T> solvable_equations =
      femib::util::remove_edges<T>(poisson.stiffness_matrix,
                                   poisson.load_vector, poisson.edge_values,
                                   not_edges);

  Eigen::Matrix<T, Eigen::Dynamic, 1> x =
      poisson.solver->solve(solvable_equations.A, solvable_equations.b);

  return femib::util::add_edges<T>(x, poisson.edge_values,
                                   poisson.V.nodes.P.size(), not_edges,
                                   poisson.V.nodes.E);
}

} // namespace femib::poisson
#endif
