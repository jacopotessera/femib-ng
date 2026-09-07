#ifndef CUDA_H_INCLUDED_
#define CUDA_H_INCLUDED_

#pragma diag_suppress = 20011, 20014, 20040

constexpr int FEMIB_CUDA_MAX_BLOCKS = 1024;
constexpr int FEMIB_CUDA_STACK_SIZE = 12928;
constexpr int FEMIB_CUDA_HEAP_SIZE = 20000000;

namespace femib::cuda {
void printSize();
size_t getStackSize();
size_t getHeapSize();
void setStackSize(size_t stackSize);
void setHeapSize(size_t heapSize);

// TODO size are probably not int
template <typename T> T *copyToDevice(T *x, int size);
template <typename T> T *copyToHost(T *X, int size);
template <typename T> void freeDevice(T *X);

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
void serial_accurate(femib::types::dvec<f, d> *X, int size_X,
                     femib::types::dtrian<f, d> *T, int size_T, bool *N);
template <typename f, int d>
void parallel_accurate(femib::types::dvec<f, d> *X, int size_X,
                       femib::types::dtrian_<f, d> *T, int size_T, bool *N);
} // namespace femib::cuda

#endif
