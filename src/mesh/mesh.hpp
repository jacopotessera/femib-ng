#ifndef MESH_HPP_INCLUDED_
#define MESH_HPP_INCLUDED_

#include "../gauss/gauss.hpp"
#include "../types/types.hpp"

namespace femib::mesh {

template <typename T, int d>
T integrate(const femib::gauss::rule<T, d> &rule,
            const std::function<T(femib::types::dvec<T, d>)> &f,
            const femib::types::dtrian<T, d> &t);
template <typename T, int d>
T integrate(const femib::gauss::rule<T, d> &rule,
            const std::function<T(femib::types::dvec<T, d>)> &f,
            const femib::types::mesh<T, d> &mesh);

} // namespace femib::mesh
#endif
