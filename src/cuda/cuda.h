#ifndef CUDA_H_INCLUDED_
#define CUDA_H_INCLUDED_

#define FEMIB_CUDA_MAX_BLOCKS 1024
#define FEMIB_CUDA_STACK_SIZE 12928
#define FEMIB_CUDA_HEAP_SIZE 20000000

namespace femib {
namespace cuda {
void printSize();
int getStackSize();
int getHeapSize();
void setStackSize(int stackSize);
void setHeapSize(int heapSize);
template <typename T> T *copyToDevice(T x);
template <typename T> T copyToHost(T *X);
} // namespace cuda
} // namespace femib

#endif
