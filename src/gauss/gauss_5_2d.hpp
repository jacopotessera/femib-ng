#ifndef GAUSS_LAGRANGE_5_2D_HPP_INCLUDED_
#define GAUSS_LAGRANGE_5_2D_HPP_INCLUDED_

#include "../types/types.hpp"

namespace femib::gauss {

// 7-point symmetric Gauss rule for the reference triangle, from Dunavant's
// tables.
template <typename T, int d> femib::gauss::rule<T, d> create_gauss_5_2d() {
  T a = 0.0597158717897698;
  T b = 0.4701420641051151;
  T wa = 0.1323941527885062 / 2.0;
  T wb = 0.1259391805448272 / 2.0;
  T wc = 0.225 / 2.0;

  femib::gauss::node<T, d> n0 = {wc, {1.0 / 3.0, 1.0 / 3.0}};
  femib::gauss::node<T, d> n1 = {wa, {a, a}};
  femib::gauss::node<T, d> n2 = {wa, {a, 1 - 2 * a}};
  femib::gauss::node<T, d> n3 = {wa, {1 - 2 * a, a}};
  femib::gauss::node<T, d> n4 = {wb, {b, b}};
  femib::gauss::node<T, d> n5 = {wb, {b, 1 - 2 * b}};
  femib::gauss::node<T, d> n6 = {wb, {1 - 2 * b, b}};

  femib::gauss::rule<T, d> rule = {{n0, n1, n2, n3, n4, n5, n6}};
  return rule;
}

} // namespace femib::gauss
#endif