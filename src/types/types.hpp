#ifndef FEMIB_TYPES_HPP
#define FEMIB_TYPES_HPP

#include <Eigen/Core>
#include <vector>

namespace femib::types {

template <typename T, int d> using dvec = Eigen::Matrix<T, d, 1>;
template <typename T, int d> using dmat = Eigen::Matrix<T, d, d>;
template <typename T, int d, int e> using rmat = Eigen::Matrix<T, d, e>;
template <typename T, int d> using dtrian = std::vector<dvec<T, d>>;
template <int d> using ditrian = Eigen::Matrix<int, d + 1, 1>;
template <int d> using ditrian_ = Eigen::Matrix<int, d, 1>;

template <typename T, int d, int e> struct xDx {
  dvec<T, e> x;
  rmat<T, d, e> dx;
};

template <typename T, int d, int e> struct F {
  std::function<dvec<T, e>(dvec<T, d>)> x;
  std::function<rmat<T, d, e>(dvec<T, d>)> dx;
  xDx<T, d, e> operator()(const dvec<T, d> &x) {
    return {this->x(x), this->dx(x)};
  };
};

template <typename f, int d> struct nodes {
  std::vector<dvec<f, d>> P;
  std::vector<std::vector<int>> T;
  std::vector<int> E;

  int get_index(int i, int n) const { return T[n][i]; }
};

template <typename f, int d> struct mesh {
  std::vector<dvec<f, d>> P;
  std::vector<ditrian<d>> T;
  std::vector<int> E;

  std::vector<dtrian<f, d>> N;
  bool initialized = false;

  void init() {
    if (!initialized) {
      for (ditrian<d> t : T) {
        dtrian<f, d> n;
        n.reserve(d + 1);
        for (int i : t) {
          n.emplace_back(P[i]);
        }
        N.emplace_back(n);
      }
    }
  }

  inline dtrian<f, d> operator[](int i) const { return N[i]; }
};

template <typename f, int d> using box = std::vector<dvec<f, d>>;

} // namespace femib::types
#endif
