#ifndef FEMIB_IB_STRUCTURE_HPP_INCLUDED_
#define FEMIB_IB_STRUCTURE_HPP_INCLUDED_

#include "../types/types.hpp"
#include <cmath>
#include <numbers>
#include <vector>

namespace femib::ib {

inline int get_next_point(const int i, const int n_points) {
  return (i + 1) % n_points;
}

template <typename T, int d> struct ring {
  std::vector<femib::types::dvec<T, d>>
      positions; // TODO we need to keep the "story"? or the story is kept in
                 // the problem?
  T spring_constant = 0;
  T reference_arc_length = 0; // circumference / n, fixed at construction
};

// TODO eh, this belongs to ring: ring.advance
template <typename T, int d>
void advance_ring(ring<T, d> &r, T time_step,
                  std::vector<femib::types::dvec<T, d>> &velocities) {
  for (size_t k = 0; k < r.positions.size(); ++k) {
    r.positions[k] += time_step * velocities[k];
  }
}

// TODO constructor?
// build a ring with center "center", radius "radius", elastic constant
// "stiffness", and "n_points" material points
template <typename T, int d>
ring<T, d> build_ring(femib::types::dvec<T, d> center, T radius, int n_points,
                      T stiffness) {
  ring<T, d> r;
  r.positions.reserve(n_points);
  const T pi = std::numbers::pi_v<T>;
  for (int i = 0; i < n_points; ++i) {
    T theta = (2 * pi * i) / static_cast<T>(n_points);
    r.positions.emplace_back(center(0) + radius * std::cos(theta),
                             center(1) + radius * std::sin(theta));
  }
  r.reference_arc_length = (2 * pi * radius) / static_cast<T>(n_points);
  // refining n_points doesn't change the physics
  r.spring_constant = stiffness / r.reference_arc_length;
  return r;
}

template <typename T, int d>
std::vector<femib::types::dvec<T, d>> elastic_force(const ring<T, d> &r) {
  int n = r.positions.size();
  std::vector<femib::types::dvec<T, d>> forces(
      n, femib::types::dvec<T, d>::Zero());
  for (int i = 0; i < n; ++i) {
    int i_next = get_next_point(i, n);
    femib::types::dvec<T, d> spring_force =
        r.spring_constant * (r.positions[i_next] - r.positions[i]);
    forces[i] += spring_force;      // pulls i toward i+1 if stretched
    forces[i_next] -= spring_force; // equal and opposite on i+1
  }
  return forces;
}

template <typename T, int d> T elastic_energy(const ring<T, d> &r) {
  const int n = r.positions.size();
  T E = 0;
  for (int i = 0; i < n; ++i) {
    int i_next = get_next_point(i, n);
    T len = (r.positions[i_next] - r.positions[i]).norm();
    E += (r.spring_constant / 2) * len * len;
  }
  return E;
}

} // namespace femib::ib
#endif
