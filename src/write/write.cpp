#include "write.hpp"

#include <array>
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
template void write_if_present<double, 1>(
    HighFive::Group &group, const std::string &name,
    const std::vector<femib::types::dvec<double, 1>> &data);
template void write_if_present<double, 2>(
    HighFive::Group &group, const std::string &name,
    const std::vector<femib::types::dvec<double, 2>> &data);

unsigned flags_for(const femib::write::mode file_mode) {
  // will fail if file already exist
  constexpr unsigned CreateNew =
      HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Excl;
  switch (file_mode) {
  case femib::write::mode::create:
    return CreateNew;
  case femib::write::mode::overwrite:
    return HighFive::File::Truncate;
  case femib::write::mode::resume:
    return HighFive::File::ReadWrite;
  }
  return CreateNew;
}
} // namespace

void femib::write::save_sim(const std::string &path,
                            const std::string &sim_name,
                            const std::string &sim_type, const mode file_mode) {
  HighFive::File file(path, flags_for(file_mode));
  if (file_mode == mode::resume) {
    return;
  }
  file.createAttribute<std::string>("sim_name", sim_name);
  file.createAttribute<std::string>("sim_type", sim_type);
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
template void femib::write::save_plot_data<double, 2>(
    const std::string &path, const femib::write::plot_data<double, 2> &data);

template <typename T>
void femib::write::save_metadata(const std::string &path,
                                 const std::string &key, const T &value) {
  HighFive::File file(path, HighFive::File::ReadWrite | HighFive::File::Create);
  HighFive::Group group = file.exist("metadata") ? file.getGroup("metadata")
                                                 : file.createGroup("metadata");
  if (group.hasAttribute(key)) {
    group.deleteAttribute(key);
  }
  group.createAttribute<T>(key, value);
}

template void femib::write::save_metadata<double>(const std::string &path,
                                                  const std::string &key,
                                                  const double &value);
template void femib::write::save_metadata<std::string>(
    const std::string &path, const std::string &key, const std::string &value);

template <typename T, int d>
void femib::write::save_mesh(const std::string &path,
                             const femib::types::mesh<T, d> &mesh) {
  HighFive::File file(path, HighFive::File::ReadWrite | HighFive::File::Create);
  HighFive::Group group = file.createGroup("mesh");

  std::vector<std::array<T, d>> points;
  points.reserve(mesh.P.size());
  for (const auto &p : mesh.P) {
    std::array<T, d> r{};
    for (int i = 0; i < d; ++i) {
      r[i] = p(i);
    }
    points.push_back(r);
  }
  group.createDataSet("points", points);

  std::vector<std::array<int, d + 1>> triangles;
  triangles.reserve(mesh.T.size());
  for (const auto &t : mesh.T) {
    std::array<int, d + 1> r{};
    for (int i = 0; i < d + 1; ++i) {
      r[i] = t(i);
    }
    triangles.push_back(r);
  }
  group.createDataSet("triangles", triangles);

  group.createDataSet("edges", mesh.E);
}

template void
femib::write::save_mesh<double, 2>(const std::string &path,
                                   const femib::types::mesh<double, 2> &mesh);
template void
femib::write::save_mesh<float, 2>(const std::string &path,
                                  const femib::types::mesh<float, 2> &mesh);
