#include "write.hpp"

#include <highfive/highfive.hpp>

namespace {

void write_if_present(HighFive::Group &group, const std::string &name,
                      const std::vector<std::vector<float>> &data) {
  if (data.empty()) {
    return;
  }
  group.createDataSet(name, data);
}

} // namespace

void femib::write::save_sim(const std::string &path,
                            const std::string &sim_name) {
  HighFive::File file(
      path, HighFive::File::Truncate); // TODO truncate? i want to continue an
                                       // existing simulation
  file.createAttribute<std::string>(
      "sim_name", sim_name); // TODO other attributes? types? finite element?
                             // mesh? delta t? etc?
}

void femib::write::save_plot_data(const std::string &path,
                                  const femib::write::plot_data &data) {
  HighFive::File file(path, HighFive::File::ReadWrite | HighFive::File::Create);
  std::string group_name = "timestep_" + std::to_string(data.time);
  HighFive::Group group = file.createGroup(group_name);

  write_if_present(group, "x", data.x);
  write_if_present(group, "u", data.u);
  write_if_present(group, "q", data.q);
}
