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

enum class mode {
  create,
  overwrite,
  resume, // TODO not yet implemented
};

void save_sim(const std::string &path, const std::string &sim_name,
              const std::string &sim_type = "", mode file_mode = mode::create);
template <typename T, int d>
void save_plot_data(const std::string &path, const plot_data<T, d> &data);

template <typename T>
void save_metadata(const std::string &path, const std::string &key,
                   const T &value);

template <typename T, int d>
void save_mesh(const std::string &path, const femib::types::mesh<T, d> &mesh);

} // namespace femib::write
#endif
