#ifndef WRITE_HPP_INCLUDED_
#define WRITE_HPP_INCLUDED_

#include <string>
#include <vector>

namespace femib::write {

struct plot_data {
  int time;
  std::vector<std::vector<float>>
      x; // N points x 2 (position) // TODO and 1? and 3?
  std::vector<std::vector<float>> u; // N points x 2 (velocity) // TODO and 3?
  std::vector<std::vector<float>> q; // M points x 1 (pressure), may be empty
};

void save_sim(const std::string &path, const std::string &sim_name);
void save_plot_data(const std::string &path, const plot_data &data);

} // namespace femib::write
#endif
