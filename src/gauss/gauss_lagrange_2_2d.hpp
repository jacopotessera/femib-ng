#ifndef GAUSS_LAGRANGE_2_2D_HPP_INCLUDED_
#define GAUSS_LAGRANGE_2_2D_HPP_INCLUDED_

#include "../types/types.hpp"

namespace femib::gauss {

template <typename T, int d>
femib::gauss::rule<T, d> create_gauss_2_2d() {
  femib::gauss::node<T, d> node1 = {1.0 / 6.0, {1.0 / 6.0, 1.0 / 6.0}};
  femib::gauss::node<T, d> node2 = {1.0 / 6.0, {1.0 / 6.0, 2.0 / 3.0}};
  femib::gauss::node<T, d> node3 = {1.0 / 6.0, {2.0 / 3.0, 1.0 / 6.0}};
  femib::gauss::rule<T, d> rule = {{node1, node2, node3}};
  return rule;
}

} // namespace femib::gauss
#endif
