#include "gauss.hpp"

#include "spdlog/spdlog.h"
#include <functional>

template <typename T, int d>
T femib::gauss::integrate(const femib::gauss::rule<T, d> rule,
                          std::function<T(femib::types::dvec<T, d>)> &f) {
  T integral = 0;
  for (femib::gauss::node<T, d> node : rule.nodes) {
    integral += node.weight * f(node.node);
    spdlog::debug("[gauss] node x  {}", node.node[0]);
    spdlog::debug("[gauss] node y  {}", node.node[1]);
    spdlog::debug("[gauss] f(node) {}", f(node.node));
  }
  return integral;
}

template float femib::gauss::integrate<float, 2>(
    const femib::gauss::rule<float, 2> rule,
    std::function<float(femib::types::dvec<float, 2>)> &f);
