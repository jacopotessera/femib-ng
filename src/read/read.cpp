#include "read.hpp"
#include "spdlog/spdlog.h"
#include <fstream>
#include <rapidcsv.h>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

template <typename T, typename W> std::vector<T> read(std::string filename);

template std::vector<femib::types::dvec<float, 2>>
read<femib::types::dvec<float, 2>, float>(std::string file);
template std::vector<femib::types::ditrian<2>>
read<femib::types::ditrian<2>, int>(std::string file);
template std::vector<int> read<int, int>(std::string file);

template <typename T, typename W> std::vector<T> read(std::string filename) {
  SPDLOG_LOGGER_DEBUG("[read] filename: {}", filename);

  rapidcsv::Document doc(filename, rapidcsv::LabelParams(-1, -1),
                         rapidcsv::SeparatorParams('\t'));

  std::vector<T> a(doc.GetRowCount());
  for (size_t i = 0; i < doc.GetRowCount(); ++i) {
    if constexpr (std::is_integral_v<T>) {
      // edges file
      a[i] = doc.GetCell<W>(0, i);
    } else {
      // points and triangles files
      // TODO replace size_t with Eigen::Index when looping on rows and cols
      for (size_t j = 0; j < a[i].size(); ++j) {
        a[i](j) = doc.GetCell<W>(j, i);
      }
    }
  }

  return a;
}

template <typename T, int d>
femib::types::mesh<T, d> read_mesh_file(std::string p, std::string t) {
  femib::types::mesh<T, d> m;
  m.P = read<femib::types::dvec<T, d>, T>(p);
  m.T = read<femib::types::ditrian<d>, int>(t);
  return m;
}

template <typename T, int d>
femib::types::mesh<T, d> read_mesh_file(std::string p, std::string t,
                                        std::string e) {
  femib::types::mesh<T, d> m = read_mesh_file<T, d>(p, t);
  m.E = read<int, int>(e);
  return m;
}

template femib::types::mesh<float, 2> read_mesh_file<float, 2>(std::string p,
                                                               std::string t);
template femib::types::mesh<float, 2>
read_mesh_file<float, 2>(std::string p, std::string t, std::string e);
