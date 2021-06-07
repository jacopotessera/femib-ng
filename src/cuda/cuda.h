#ifndef CUDA_H_INCLUDED_
#define CUDA_H_INCLUDED_

constexpr int FEMIB_CUDA_MAX_BLOCKS = 1024;
constexpr int FEMIB_CUDA_STACK_SIZE = 12928;
constexpr int FEMIB_CUDA_HEAP_SIZE = 20000000;

namespace femib::cuda {
void printSize();
int getStackSize();
int getHeapSize();
void setStackSize(int stackSize);
void setHeapSize(int heapSize);
template <typename T> T *copyToDevice(T x);
template <typename T> T copyToHost(T *X);

template <typename f, int d>
bool in_box(const femib::types::dvec<f, d> &P,
            const femib::types::dtrian<f, d> &T);
template <typename f, int d>
bool in_triangle(const femib::types::dvec<f, d> &P,
                 const femib::types::dtrian<f, d> &T);

} // namespace femib::cuda

#endif
