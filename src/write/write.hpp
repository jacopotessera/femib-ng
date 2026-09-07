#ifndef WRITE_HPP_INCLUDED_
#define WRITE_HPP_INCLUDED_

#include <string>
#include <vector>

#include "../types/types.hpp"

namespace femib::write {

template <typename T, int d> struct plot_data {
  int time;
  std::vector<femib::types::dvec<T, d>> x; // position
  std::vector<femib::types::dvec<T, d>> u; // velocity
  std::vector<femib::types::dvec<T, 1>> q; // pressure
  std::vector<femib::types::dvec<T, d>> X; // structure
};

void save_sim(const std::string &path, const std::string &sim_name);
template <typename T, int d>
void save_plot_data(const std::string &path, const plot_data<T, d> &data);

} // namespace femib::write
#endif
