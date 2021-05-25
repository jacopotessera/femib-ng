#ifndef GAUSS_HPP_INCLUDED_
#define GAUSS_HPP_INCLUDED_

#include "../types/types.hpp"
#include <execution>
#include <functional>
#include <numeric>
#include <vector>

namespace femib::gauss {

template <typename T, int d> struct node {
  T weight;
  femib::types::dvec<T, d> node;
};

template <typename T, int d> struct rule { std::vector<node<T, d>> nodes; };

template <typename T, int d>
T integrate(const femib::gauss::rule<T, d> rule,
            std::function<T(femib::types::dvec<T, d>)> &f) {
  auto unary_op = [&f](const femib::gauss::node<T, d> &node) {
    return node.weight * f(node.node);
  };
  return std::transform_reduce(std::execution::seq, rule.nodes.begin(),
                               rule.nodes.end(), 0.0, std::plus<>(), unary_op);
}

} // namespace femib::gauss
#endif
