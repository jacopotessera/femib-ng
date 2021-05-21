#ifndef FEMIB_TYPES_HPP
#define FEMIB_TYPES_HPP

#include "iterator_tpl.h"
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
};

/// template <typename T, int d> using mesh = std::vector<dtrian<T, d>>;
template <typename f, int d> struct mesh {
  std::vector<dvec<f, d>> P;
  std::vector<ditrian<d>> T;
  std::vector<int> E;

  typedef dtrian<f, d> dtrianfd;
  inline dtrianfd operator[](int pos) {
    dtrianfd t;
    t.reserve(d + 1);
    for (int tt : T[pos]) {
      t.emplace_back(P[tt]);
    }
    return t;
  }

  typedef ditrian<d> ditriand;
  STL_TYPEDEFS(ditriand);
  struct it_state {
    int pos;
    inline void next(const mesh<f, d> *ref) { ++pos; }
    inline void begin(const mesh<f, d> *ref) { pos = 0; }
    inline void end(const mesh<f, d> *ref) { pos = ref->T.size(); }
    inline dtrianfd get(mesh<f, d> *ref) {
      dtrianfd t;
      t.reserve(d + 1);
      for (int tt : ref->T[pos]) {
        t.emplace_back(ref->P[tt]);
      }
      return t;
      // return
      // {ref->P[(ref->T[pos])[0]],ref->P[(ref->T[pos])[1]],ref->P[(ref->T[pos])[2]]};
    }
    inline dtrianfd get(const mesh<f, d> *ref) {
      dtrianfd t;
      t.reserve(d + 1);
      for (int tt : ref->T[pos]) {
        t.emplace_back(ref->P[tt]);
      }
      return t;
      // return
      // {ref->P[(ref->T[pos])[0]],ref->P[(ref->T[pos])[1]],ref->P[(ref->T[pos])[2]]};
    }
    inline bool cmp(const it_state &s) const { return pos != s.pos; }
  };
  SETUP_ITERATORS(mesh, dtrianfd, it_state);
};
} // namespace femib::types
#endif
