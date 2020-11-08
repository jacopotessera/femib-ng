#ifndef FEMIB_TYPES_HPP
#define FEMIB_TYPES_HPP

#include "iterator_tpl.h"
#include <Eigen/Core>
#include <vector>

namespace femib::types {

template <typename T, int d> using dvec = Eigen::Matrix<T, d, 1>;
template <typename T, int d> using dmat = Eigen::Matrix<T, d, d>;
template <typename T, int d> using dtrian = std::vector<dvec<T, d>>;
template <int d> using ditrian = Eigen::Matrix<int, d + 1, 1>;

template <typename T, int d> struct xDx {
  dvec<T, d> x;
  dmat<T, d> dx;
};

template <typename T, int d> struct F {
  std::function<dvec<T, d>(dvec<T, d>)> x;
  std::function<dmat<T, d>(dvec<T, d>)> dx;
  xDx<T, d> operator()(const dvec<T, d> &x) {
    return {this->x(x), this->dx(x)};
  };
};

/// template <typename T, int d> using mesh = std::vector<dtrian<T, d>>;
template <typename f, int d> struct mesh {
  std::vector<dvec<f, d>> P;
  std::vector<ditrian<d>> T;
  std::vector<int> E;

  typedef ditrian<d> ditriand;
  typedef dtrian<f, d> dtrianfd;
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
    inline const dtrianfd get(const mesh<f, d> *ref) {
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
