#ifndef CUDA_H_INCLUDED_
#define CUDA_H_INCLUDED_

#include "../types/types.hpp"

#pragma diag_suppress = 20011, 20014, 20040

constexpr int FEMIB_CUDA_MAX_BLOCKS = 1024;
constexpr int FEMIB_CUDA_STACK_SIZE = 12928;
constexpr int FEMIB_CUDA_HEAP_SIZE = 20000000;

namespace femib::cuda {
void print_size();
size_t get_stack_size();
size_t get_heap_size();
void set_stack_size(size_t stack_size);
void set_heap_size(size_t heap_size);

template <typename T> T *copy_to_device(T *x, size_t size);
template <typename T> T *copy_to_host(T *X, size_t size);
template <typename T> void free_device(T *X);
template <typename T> T *alloc_device(size_t size);
template <typename T> void free_device_quiet(T *X);

template <typename f, int d>
bool in_box(const femib::types::dvec<f, d> &P,
            const femib::types::dtrian<f, d> &T);
template <typename f, int d>
bool in_triangle(const femib::types::dvec<f, d> &P,
                 const femib::types::dtrian<f, d> &T);
template <typename f, int d>
bool accurate(const femib::types::dvec<f, d> &P,
              const femib::types::dtrian<f, d> &T);
template <typename f, int d>
void serial_accurate(femib::types::dvec<f, d> *X, size_t size_X,
                     femib::types::dtrian<f, d> *T, size_t size_T, bool *N);
// throws std::invalid_argument if size_X or size_T exceed the device's grid
// limits
template <typename f, int d>
void parallel_accurate(femib::types::dvec<f, d> *X, size_t size_X,
                       femib::types::dtrian_<f, d> *T, size_t size_T, bool *N);
} // namespace femib::cuda

#endif
