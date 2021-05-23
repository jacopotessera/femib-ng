#include "gauss.hpp"

#include "spdlog/spdlog.h"
#include <execution>
#include <functional>
#include <numeric>

template <typename T, int d>
T femib::gauss::integrate(const femib::gauss::rule<T, d> rule,
                          std::function<T(femib::types::dvec<T, d>)> &f) {
  auto unary_op = [&f](const femib::gauss::node<T, d> &node) {
    return node.weight * f(node.node);
  };
  return std::transform_reduce(std::execution::seq, rule.nodes.begin(),
                               rule.nodes.end(), 0.0, std::plus<>(), unary_op);
}

template float femib::gauss::integrate<float, 2>(
    const femib::gauss::rule<float, 2> rule,
    std::function<float(femib::types::dvec<float, 2>)> &f);
