#ifndef FEMIB_HPP_INCLUDED_
#define FEMIB_HPP_INCLUDED_

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <execution>
#include <numeric>
#include <vector>

#include "../affine/affine.hpp"
#include "../finite_element_space/finite_element_space.hpp"
#include "../gauss/gauss.hpp"
#include "../mesh/mesh.hpp"
#include "../types/differential_operation.hpp"
#include "../types/types.hpp"

namespace femib::util {

template <typename T>
Eigen::SparseMatrix<T>
triplets2sparse(const std::vector<Eigen::Triplet<T>> &triplets, const int rows,
                const int cols) {
  Eigen::SparseMatrix<T> sparse_matrix = Eigen::SparseMatrix<T>(rows, cols);
  sparse_matrix.setFromTriplets(triplets.begin(), triplets.end());
  return sparse_matrix;
}
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
triplets2dense(const std::vector<Eigen::Triplet<T>> &triplets, const int rows,
               const int cols) {
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> dense_matrix =
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(
          triplets2sparse<T>(triplets, rows, cols));
  return dense_matrix;
}

// selects (and permutes) the given rows/cols out of a matrix: P * A * P^T
// keeps only rows/cols in `rows`, in that order
template <typename T>
Eigen::SparseMatrix<T> selection_matrix(const std::vector<int> &rows, int n) {
  Eigen::SparseMatrix<T> P(rows.size(), n);
  std::vector<Eigen::Triplet<T>> triplets;
  triplets.reserve(rows.size());
  for (size_t i = 0; i < rows.size(); ++i) {
    triplets.push_back(Eigen::Triplet<T>((int)i, rows[i], T(1)));
  }
  P.setFromTriplets(triplets.begin(), triplets.end());
  return P;
}

// TODO are the captures correct?
template <typename T, int d, int e>
femib::types::F<T, d, e> base_function2real_function(
    const femib::finite_element_space::finite_element_space<T, d, e> &v, int i,
    const femib::types::dmat<T, d> &Binv, const femib::types::dvec<T, d> &b) {
  femib::types::F<T, d, e> a;
  a.x =
      [&v, i, Binv, b](const femib::types::dvec<T, d> &x) {
        return (v.finite_element.base_functions[i].x(Binv * (x - b)));
      },
  a.dx = [&v, i, Binv,
          b](const femib::types::dvec<T, d> &x) -> femib::types::rmat<T, d, e> {
    return (Binv.transpose() *
            v.finite_element.base_functions[i].dx(Binv * (x - b)));
  };
  return a;
}

template <typename T> struct sparse_solvable_equations {
  Eigen::SparseMatrix<T> A;
  Eigen::Matrix<T, Eigen::Dynamic, 1> b;
};

// reduces A x = f (with Dirichlet values bV on the E nodes) to a smaller
// system over just not_edges, folding bV's contribution into the RHS.
template <typename T>
sparse_solvable_equations<T>
remove_edges(const Eigen::SparseMatrix<T> &dM,              // TODO full matrix
             const Eigen::Matrix<T, Eigen::Dynamic, 1> &dF, // TODO full vector
             const Eigen::Matrix<T, Eigen::Dynamic, 1> &bV, // TODO edge values?
             const std::vector<int> &not_edges) { // TODO not edges nodes

  Eigen::Matrix<T, Eigen::Dynamic, 1> ss = dM * bV;

  int n = static_cast<int>(not_edges.size());
  Eigen::Matrix<T, Eigen::Dynamic, 1> bbb(n);
  for (int i = 0; i < n; ++i) {
    bbb(i) = dF(not_edges[i]) - ss(not_edges[i]);
  }

  Eigen::SparseMatrix<T> P =
      selection_matrix<T>(not_edges, static_cast<int>(dM.rows()));
  Eigen::SparseMatrix<T> AAA = P * dM * P.transpose();

  return {AAA, bbb};
}

// inverse of remove_edges: expands a not_edges-only solution xxx back to the
// full `rows`-length vector, filling E nodes with bV and everything else
// with xxx's corresponding not_edges entry.
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1>
add_edges(const Eigen::Matrix<T, Eigen::Dynamic, 1>
              &xxx, // TODO solution for the reduced problem
          const Eigen::Matrix<T, Eigen::Dynamic, 1> &bV,
          int rows, // TODO edge values
          const std::vector<int> &not_edges, const std::vector<int> &nodesE) {

  Eigen::Matrix<T, Eigen::Dynamic, 1> xx; // TODO full solution vector
  xx.resize(rows, 1);

  for (int i = 0; i < rows; i++) {
    xx(i, 0) = 0.0;
    auto k = std::find(not_edges.begin(), not_edges.end(), i);
    if (k != not_edges.end()) {
      xx(i, 0) = xxx(k - not_edges.begin(), 0);
    }
    auto kk = std::find(nodesE.begin(), nodesE.end(), i);
    if (kk != nodesE.end()) {
      xx(i, 0) = bV(i, 0);
    }
  }

  return xx;
}

template <typename T, int d, int e>
std::vector<Eigen::Triplet<T>> build_diagonal_matrix(
    const femib::finite_element_space::finite_element_space<T, d, e> &v,
    const femib::gauss::rule<T, d> &integration_rule,
    const std::function<std::function<T(femib::types::dvec<T, d>)>(
        femib::types::F<T, d, e>, femib::types::F<T, d, e>)> &_operator) {
  const int basis_count =
      static_cast<int>(v.finite_element.base_functions.size());
  const int n_tri = static_cast<int>(v.mesh.T.size());

  std::vector<Eigen::Triplet<T>> diagonal_matrix(static_cast<size_t>(n_tri) *
                                                 basis_count * basis_count);

  std::vector<int> tri_indices(n_tri);
  std::iota(tri_indices.begin(), tri_indices.end(), 0);

  std::for_each(
      std::execution::par, tri_indices.begin(), tri_indices.end(), [&](int n) {
        femib::types::dtrian<T, d> triangle = v.mesh[n];
        femib::types::dmat<T, d> affine_Binv =
            femib::affine::affineBinv(triangle);
        femib::types::dvec<T, d> affine_b = femib::affine::affineb(triangle);
        for (int i = 0; i < basis_count; ++i) {
          femib::types::F<T, d, e> real_function_i =
              femib::util::base_function2real_function<T, d, e>(
                  v, i, affine_Binv, affine_b);
          for (int j = 0; j < basis_count; ++j) {
            femib::types::F<T, d, e> real_function_j =
                femib::util::base_function2real_function<T, d, e>(
                    v, j, affine_Binv, affine_b);
            T val = femib::mesh::integrate<T, d>(
                integration_rule, _operator(real_function_i, real_function_j),
                triangle);
            size_t idx =
                (static_cast<size_t>(n) * basis_count + i) * basis_count + j;
            diagonal_matrix[idx] = Eigen::Triplet<T>(
                v.nodes.get_index(i, n), v.nodes.get_index(j, n),
                val); // TODO thread-safe? yes, but...
          }
        }
      });
  return diagonal_matrix;
}

template <typename T, int d, int e>
Eigen::SparseMatrix<T> build_stiffness_matrix(
    const femib::finite_element_space::finite_element_space<T, d, e> &v,
    const femib::gauss::rule<T, d> &integration_rule) {
  return triplets2sparse(
      build_diagonal_matrix(v, integration_rule, ddot<T, d, e>),
      v.nodes.P.size(), v.nodes.P.size());
}

template <typename T, int d, int e>
Eigen::SparseMatrix<T> build_mass_matrix(
    const femib::finite_element_space::finite_element_space<T, d, e> &v,
    const femib::gauss::rule<T, d> &integration_rule) {
  return triplets2sparse(
      build_diagonal_matrix(v, integration_rule, mass<T, d, e>),
      v.nodes.P.size(),
      v.nodes.P.size()); // TODO add v.size()
}

template <typename T, int d, int e>
Eigen::Matrix<T, Eigen::Dynamic, 1> build_load_vector(
    const femib::finite_element_space::finite_element_space<T, d, e> &v,
    const femib::gauss::rule<T, d> &quadrature_rule,
    const std::function<std::function<T(femib::types::dvec<T, d>)>(
        femib::types::F<T, d, e>)> &body_force) { // TODO body force?
  const int basis_count =
      static_cast<int>(v.finite_element.base_functions.size());
  const int n_tri = static_cast<int>(v.mesh.T.size());

  std::vector<Eigen::Triplet<T>> load_vector(static_cast<size_t>(n_tri) *
                                             basis_count);

  std::vector<int> tri_indices(n_tri);
  std::iota(tri_indices.begin(), tri_indices.end(), 0); // 0 -> n_tri - 1

  std::for_each(
      std::execution::par, tri_indices.begin(), tri_indices.end(), [&](int n) {
        femib::types::dtrian<T, d> triangle = v.mesh[n];
        femib::types::dmat<T, d> affine_Binv =
            femib::affine::affineBinv(triangle);
        femib::types::dvec<T, d> affine_b = femib::affine::affineb(triangle);
        for (int i = 0; i < basis_count; ++i) {
          femib::types::F<T, d, e> real_function =
              femib::util::base_function2real_function<T, d, e>(
                  v, i, affine_Binv, affine_b);
          T val = femib::mesh::integrate<T, d>(
              quadrature_rule, body_force(real_function), triangle);
          size_t idx = static_cast<size_t>(n) * basis_count + i;
          load_vector[idx] = Eigen::Triplet<T>(
              v.nodes.get_index(i, n), 0, val); // TODO thread-safe? yes, but...
        }
      });
  return Eigen::Matrix<T, Eigen::Dynamic, 1>(
      triplets2dense(load_vector, v.nodes.P.size(), 1));
}

template <typename T, int d>
std::vector<Eigen::Triplet<T>> build_off_diagonal_matrix(
    const femib::finite_element_space::finite_element_space<T, d, d> &v,
    const femib::finite_element_space::finite_element_space<T, d, 1> &q,
    const femib::gauss::rule<T, d> &rule,
    const std::function<std::function<T(femib::types::dvec<T, d>)>(
        femib::types::F<T, d, d>, femib::types::F<T, d, 1>)> &fff) {
  std::vector<Eigen::Triplet<T>> BB;
  for (int n = 0; n < v.mesh.T.size(); ++n) {
    femib::types::dtrian<T, d> t = v.mesh[n];
    femib::types::dmat<T, d> Binv = femib::affine::affineBinv(t);
    femib::types::dvec<T, d> bb = femib::affine::affineb(t);
    for (int i = 0; i < v.finite_element.base_functions.size(); ++i) {
      femib::types::F<T, d, d> a =
          femib::util::base_function2real_function<T, d, d>(v, i, Binv, bb);
      for (int j = 0; j < q.finite_element.base_functions.size(); ++j) {
        femib::types::F<T, d, 1> b =
            femib::util::base_function2real_function<T, d, 1>(q, j, Binv, bb);
        T m = femib::mesh::integrate<T, d>(rule, fff(a, b), t);
        BB.push_back(Eigen::Triplet<T>(v.nodes.get_index(i, n),
                                       q.nodes.get_index(j, n), m));
      }
    }
  }
  return BB;
}

template <typename T, int d, int e>
std::vector<Eigen::Triplet<T>>
build_edges(const femib::finite_element_space::finite_element_space<T, d, e> &s,
            const std::function<T(femib::types::dvec<T, d>)> &b) {
  std::vector<Eigen::Triplet<T>> B;
  for (int i : s.nodes.E) {
    B.push_back(Eigen::Triplet<T>(i, 0, b(s.nodes.P[i])));
  }
  return B;
}

// TODO un-normalized pressure constraint
// TODO what´s the difference with build_load_vector
template <typename T, int d>
Eigen::Matrix<T, 1, Eigen::Dynamic> build_domain_integral_row(
    const femib::finite_element_space::finite_element_space<T, d, 1> &v,
    const femib::gauss::rule<T, d> &rule) {
  std::vector<Eigen::Triplet<T>> B;
  for (int n = 0; n < v.mesh.T.size(); ++n) {
    femib::types::dtrian<T, d> t = v.mesh[n];
    femib::types::dmat<T, d> Binv = femib::affine::affineBinv(t);
    femib::types::dvec<T, d> bb = femib::affine::affineb(t);
    for (int i = 0; i < v.finite_element.base_functions.size(); ++i) {
      femib::types::F<T, d, 1> a =
          femib::util::base_function2real_function<T, d, 1>(v, i, Binv, bb);
      auto g = [&](const femib::types::dvec<T, d> &x) { return a.x(x)(0); };
      T m = femib::mesh::integrate<T, d>(rule, g, t);
      B.push_back(Eigen::Triplet<T>(0, v.nodes.get_index(i, n), m));
    }
  }
  return femib::util::triplets2dense<T>(B, 1, v.nodes.P.size());
}

template <typename T, int d, int e>
std::vector<int> build_not_edges(
    const femib::finite_element_space::finite_element_space<T, d, e> &s) {
  std::vector<bool> is_edge(s.nodes.P.size(), false);
  for (int i : s.nodes.E) {
    is_edge[i] = true;
  }
  std::vector<int> not_edges;
  not_edges.reserve(s.nodes.P.size());
  for (int i = 0; i < static_cast<int>(s.nodes.P.size()); i++) {
    if (!is_edge[i]) {
      not_edges.push_back(i);
    }
  }
  return not_edges;
}

} // namespace femib::util
#endif
