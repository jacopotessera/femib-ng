#ifndef FEMIB_IB_STRUCTURE_HPP_INCLUDED_
#define FEMIB_IB_STRUCTURE_HPP_INCLUDED_

#include "../finite_element_space/finite_element_space.hpp"
#include "../types/types.hpp"
#include <cmath>
#include <numbers>
#include <vector>

namespace femib::ib {

inline int get_next_point(const int i, const int n_points) {
  return (i + 1) % n_points;
}

template <typename T, int d> struct ring {
  std::vector<femib::types::dvec<T, d>> X; // TODO we need to keep the "story"
  T k = 0;
  T dS = 0; // fixed material arc-length weight per point (circumference / n)
};

// TODO eh, this belongs to ring: ring.advance
template <typename T, int d>
void advance_ring(ring<T, d> &r, T deltat,
                  std::vector<femib::types::dvec<T, d>> &U) {
  for (size_t k = 0; k < r.X.size(); ++k) {
    r.X[k] += deltat * U[k];
  }
}

// build a ring with center "c", radius "radius", elastic constant "k", and
// "n_points" material points
template <typename T, int d>
ring<T, d> build_ring(femib::types::dvec<T, d> c, T radius, int n_points, T k) {
  ring<T, d> r;
  r.k = k;
  r.X.reserve(n_points);
  const T pi = std::numbers::pi_v<T>;
  for (int i = 0; i < n_points; ++i) {
    T theta = (2 * pi * i) / static_cast<T>(n_points);
    r.X.emplace_back(c(0) + radius * std::cos(theta),
                     c(1) + radius * std::sin(theta));
  }
  r.dS = (2 * pi * radius) / static_cast<T>(n_points);
  return r;
}

template <typename T, int d>
std::vector<femib::types::dvec<T, d>> elastic_force(const ring<T, d> &r) {
  int n = r.X.size();
  std::vector<femib::types::dvec<T, d>> F(n, femib::types::dvec<T, d>::Zero());
  for (int i = 0; i < n; ++i) {
    int i_next = get_next_point(i, n);
    femib::types::dvec<T, d> spring_force = r.k * (r.X[i_next] - r.X[i]);
    F[i] += spring_force;      // pulls i toward i+1 if stretched
    F[i_next] -= spring_force; // equal and opposite on i+1
  }
  return F;
}

template <typename T, int d>
Eigen::Matrix<T, Eigen::Dynamic, 1>
spread_force(ring<T, d> &r,
             femib::finite_element_space::finite_element_space<T, d, d> &V) {
  std::vector<types::dvec<T, d>> forces = elastic_force(r);
  return spread_force<T, d>(V, r.X, forces, r.dS);
}

template <typename T, int d> T elastic_energy(const ring<T, d> &r) {
  const int n = r.X.size();
  T E = 0;
  for (int i = 0; i < n; ++i) {
    int i_next = get_next_point(i, n);
    T len = (r.X[i_next] - r.X[i]).norm();
    E += (r.k / 2) * len * len;
  }
  return E;
}

} // namespace femib::ib
#endif
