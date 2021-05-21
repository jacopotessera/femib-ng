/*
 *	read.h
 */

#ifndef READ_H_INCLUDED_
#define READ_H_INCLUDED_

#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

template <class T, class W> std::vector<T> read(std::string file);

template <class W> W castToT(std::string s);

template <class T, class W> void set(T &a, int i, std::string s);

femib::types::mesh<float, 2> readMesh(std::string p, std::string t);
femib::types::mesh<float, 2> readMesh(std::string p, std::string t,
                                      std::string e);

#endif
