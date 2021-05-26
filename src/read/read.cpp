#include "read.hpp"
#include "spdlog/spdlog.h"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

template <class T, class W> void set(T &a, int i, std::string token);
template <class T> T castToT(std::string s);
template <class W> void castToT(std::string s, W &w);
template <typename T, typename W> std::vector<T> read(std::string filename);

// template <> double castToT<double>(std::string s) { return std::stod(s); }
template <> float castToT<float>(std::string s) { return std::stof(s); }
template <> int castToT<int>(std::string s) { return std::stoi(s); }

template <class T, class W> void set(T &a, int i, std::string token) {
  a(i) = castToT<W>(token);
}

template <> void set<int, int>(int &a, int i, std::string token) {
  a = castToT<int>(token);
}

template std::vector<femib::types::dvec<float, 2>>
read<femib::types::dvec<float, 2>, float>(std::string file);
template std::vector<femib::types::ditrian<2>>
read<femib::types::ditrian<2>, int>(std::string file);
template std::vector<int> read<int, int>(std::string file);

// template double castToT<double>(std::string file);
template float castToT<float>(std::string file);
template int castToT<int>(std::string file);

template void
set<femib::types::dvec<float, 2>, float>(femib::types::dvec<float, 2> &a, int i,
                                         std::string token);
template void set<femib::types::ditrian<2>, int>(femib::types::ditrian<2> &a,
                                                 int i, std::string token);
template void set<int, int>(int &a, int i, std::string token);


template <typename T, typename W> std::vector<T> read(std::string filename) {
  std::string tab = "\t";
  std::vector<T> a;
  std::string line;
  std::ifstream file(filename);

  spdlog::debug("[read] filename: {}", filename);
  if (file.is_open()) {
    for (int i = 0; std::getline(file, line); ++i) {
      // spdlog::debug("[read] line  {}", i);
      T t;
      a.push_back(t);
      size_t pos = 0;
      std::string token;
      for (int j = 0; (pos = line.find(tab)) != std::string::npos; ++j) {
        token = line.substr(0, pos);
        line.erase(0, pos + tab.length());
        // spdlog::debug("[read] token  {}", token);
        set<T, W>(a[i], j, token);
      }
    }
    file.close();
  } else {
    throw std::runtime_error("Unable to open file " + filename);
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
  femib::types::mesh<T, d> m = read_mesh_file<T,d>(p, t);
  m.E = read<int, int>(e);
  return m;
}

template
femib::types::mesh<float, 2> read_mesh_file<float, 2>(std::string p, std::string t);
template
femib::types::mesh<float, 2> read_mesh_file<float, 2>(std::string p, std::string t, std::string e);
