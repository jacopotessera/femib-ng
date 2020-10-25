#ifndef GAUSS_HPP_INCLUDED_
#define GAUSS_HPP_INCLUDED_

#include "../types/types.hpp"
#include <functional>
#include <vector>

namespace femib::gauss {

template <typename T, int d> struct node {
  T weight;
  femib::types::dvec<T, d> node;
};

template <typename T, int d> struct rule { std::vector<node<T, d>> nodes; };

template <typename T, int d>
T integrate(const rule<T, d> rule,
            std::function<T(femib::types::dvec<T, d>)> &f);

} // namespace femib::gauss
#endif
