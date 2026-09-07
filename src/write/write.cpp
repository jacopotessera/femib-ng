#include "write.hpp"

#include <highfive/highfive.hpp>

namespace {
template <typename T, int d>
void write_if_present(HighFive::Group &group, const std::string &name,
                      const std::vector<femib::types::dvec<T, d>> &data) {
  if (data.empty()) {
    return;
  }
  std::vector<std::array<T, d>> records;
  records.reserve(data.size());
  for (const auto &v : data) {
    std::array<T, d> r{};
    for (int i = 0; i < d; ++i) {
      r[i] = v(i);
    }
    records.push_back(r);
  }
  group.createDataSet(name, records);
}

template void write_if_present<float, 1>(
    HighFive::Group &group, const std::string &name,
    const std::vector<femib::types::dvec<float, 1>> &data);
template void write_if_present<float, 2>(
    HighFive::Group &group, const std::string &name,
    const std::vector<femib::types::dvec<float, 2>> &data);
} // namespace

void femib::write::save_sim(const std::string &path,
                            const std::string &sim_name) {
  HighFive::File file(
      path, HighFive::File::Truncate); // TODO truncate? i want to continue an
                                       //  existing simulation
  file.createAttribute<std::string>(
      "sim_name", sim_name); // TODO other attributes? types? finite element?
                             //  mesh? delta t? etc?
}

template <typename T, int d>
void femib::write::save_plot_data(const std::string &path,
                                  const femib::write::plot_data<T, d> &data) {
  HighFive::File file(path, HighFive::File::ReadWrite | HighFive::File::Create);
  const std::string group_name = "timestep_" + std::to_string(data.time);
  HighFive::Group group = file.createGroup(group_name);

  write_if_present(group, "x", data.x);
  write_if_present(group, "u", data.u);
  write_if_present(group, "q", data.q);
  write_if_present(group, "X", data.X);
}

template void femib::write::save_plot_data<float, 2>(
    const std::string &path, const femib::write::plot_data<float, 2> &data);