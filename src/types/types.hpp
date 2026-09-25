#ifndef FEMIB_TYPES_HPP
#define FEMIB_TYPES_HPP

#include "spdlog/spdlog.h"
#include <Eigen/Core>
#include <memory>
#include <vector>

namespace femib::types {

// TODO we like these names?
template <typename T, int d> using dvec = Eigen::Matrix<T, d, 1>;
template <typename T, int d> using dmat = Eigen::Matrix<T, d, d>;
template <typename T, int d, int e> using rmat = Eigen::Matrix<T, d, e>;
template <typename T, int d> using dtrian = std::vector<dvec<T, d>>;
template <typename T, int d> using dtrian_ = dvec<T, d>[d + 1];
template <int d> using ditrian = Eigen::Matrix<int, d + 1, 1>;
template <int d> using ditrian_ = Eigen::Matrix<int, d, 1>;

template <typename T, int d, int e> struct xDx {
  dvec<T, e> x;
  rmat<T, d, e> dx;
};

template <typename T, int d, int e> struct F {
  std::function<dvec<T, e>(dvec<T, d>)> x;
  std::function<rmat<T, d, e>(dvec<T, d>)> dx;
  xDx<T, d, e> operator()(const dvec<T, d> &p) {
    return {this->x(p), this->dx(p)};
  };
};

template <typename f, int d> struct nodes {
  std::vector<dvec<f, d>> P;
  std::vector<std::vector<int>> T;
  std::vector<int> E;

  [[nodiscard]] int get_index(const int i, const int n) const {
    return T[n][i];
  }
};

template <typename f, int d> struct mesh {
  std::vector<dvec<f, d>> P;
  std::vector<ditrian<d>> T;
  std::vector<int> E;

  std::vector<dtrian<f, d>> N;

  mutable std::shared_ptr<void> device_triangles_cache;

  mesh() = default;
  mesh(std::vector<dvec<f, d>> p, std::vector<ditrian<d>> t,
       std::vector<int> e = {})
      : P(std::move(p)), T(std::move(t)), E(std::move(e)) {
    N.reserve(T.size());
    for (const ditrian<d> &t_ : T) {
      dtrian<f, d> n;
      n.reserve(d + 1);
      for (int i : t_) {
        n.emplace_back(P[i]);
      }
      N.emplace_back(std::move(n));
    }
  }

  inline const dtrian<f, d> &operator[](int i) const { return N[i]; }
};

template <typename f, int d> struct box {
  dvec<f, d> bottom;
  dvec<f, d> top;
};

template <typename f, int d>
dtrian_<f, d> *
vector_dtrian2pointer_dtrian_(const std::vector<dtrian<f, d>> &N) {
  auto *p = new dtrian_<f, d>[N.size()];
  for (int i = 0; i < N.size(); ++i) {
    p[i][0] = N[i][0];
    p[i][1] = N[i][1];
    p[i][2] = N[i][2];
  }
  return p;
}

} // namespace femib::types
#endif
