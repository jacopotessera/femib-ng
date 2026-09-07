#ifndef FEMIB_UTILS_HPP
#define FEMIB_UTILS_HPP

#include "types.hpp"

inline femib::types::mesh<float, 2> make_unit_square_mesh(const int n) {
  femib::types::mesh<float, 2> mesh;
  int side = n + 1;
  auto idx = [side](const int i, const int j) { return i * side + j; };
  for (int i = 0; i < side; ++i)
    for (int j = 0; j < side; ++j) {
      float x = static_cast<float>(i) / static_cast<float>(n);
      float y = static_cast<float>(j) / static_cast<float>(n);
      mesh.P.emplace_back(x, y);
    }
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) {
      int p00 = idx(i, j), p10 = idx(i + 1, j), p01 = idx(i, j + 1),
          p11 = idx(i + 1, j + 1);
      mesh.T.emplace_back(p00, p10, p11);
      mesh.T.emplace_back(p00, p11, p01);
    }
  for (int i = 0; i < side; ++i)
    for (int j = 0; j < side; ++j)
      if (i == 0 || i == n || j == 0 || j == n)
        mesh.E.push_back(idx(i, j));
  return mesh;
}

inline std::string get_time() {
  const auto now = std::chrono::system_clock::now();
  return std::format("{0:%F}T{0:%T}", now);
}

#endif // FEMIB_UTILS_HPP
