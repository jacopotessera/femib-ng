#ifndef FINITE_ELEMENT_SPACE_HPP_INCLUDED_
#define FINITE_ELEMENT_SPACE_HPP_INCLUDED_
#include "../types/types.hpp"
#include "../finite_element/finite_element.hpp"

namespace femib::finite_element_space {

  template <typename T, int d, int e>
  struct finite_element_space {
	  femib::finite_element::finite_element<T,d,e> finite_element;
	  femib::types::mesh<T,d> mesh;
  };

}
#endif
