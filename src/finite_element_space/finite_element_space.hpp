#ifndef FINITE_ELEMENT_SPACE_HPP_INCLUDED_
#define FINITE_ELEMENT_SPACE_HPP_INCLUDED_
#include "../finite_element/finite_element.hpp"
#include "../types/types.hpp"

namespace femib::finite_element_space {

template <typename T, int d, int e> struct finite_element_space {
  femib::finite_element::finite_element<T, d, e> finite_element;
  femib::types::mesh<T, d> mesh;
  femib::types::nodes<T, d> nodes;
};

} // namespace femib::finite_element_space
#endif
