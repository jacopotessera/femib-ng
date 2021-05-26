/*
 *	read.h
 */

#ifndef READ_H_INCLUDED_
#define READ_H_INCLUDED_

#include "../types/types.hpp"

// namespace femib::read {

template <typename T, int d>
femib::types::mesh<T, d> read_mesh_file(std::string p, std::string t);

template <typename T, int d>
femib::types::mesh<T, d> read_mesh_file(std::string p, std::string t,
                                        std::string e);
//}
#endif
