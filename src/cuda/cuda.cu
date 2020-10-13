#include <iostream>

#include "cuda.h"
#include "spdlog/spdlog.h"
#include "mini-book.h"


void femib::cuda::printSize() {
	spdlog::set_pattern("[%Y-%m-%dT%T] [%l] [%@@%!] %v");
	SPDLOG_INFO("[CUDA stack size] found to be {}", getStackSize());
	SPDLOG_INFO("[CUDA heap  size] found to be {}", getHeapSize());
}

int femib::cuda::getStackSize() {
	size_t size_stack;
	cudaDeviceGetLimit(&size_stack, cudaLimitStackSize);
	return (int)size_stack;
}

int femib::cuda::getHeapSize() {
	size_t size_heap;
	cudaDeviceGetLimit(&size_heap, cudaLimitMallocHeapSize);
	return (int)size_heap;
}

void femib::cuda::setStackSize(int stackSize) {
	cudaDeviceSetLimit(cudaLimitStackSize,stackSize);
}

void femib::cuda::setHeapSize(int heapSize) {
	cudaDeviceSetLimit(cudaLimitMallocHeapSize,heapSize*sizeof(double));
}

template<typename T>
T* femib::cuda::copyToDevice(T x) {
	T *X;
    HANDLE_ERROR(cudaMalloc((void**)&X,sizeof(T)));
    HANDLE_ERROR(cudaMemcpy(X,&x,sizeof(T),cudaMemcpyHostToDevice));
    return X;
}

template<typename T>
T femib::cuda::copyToHost(T *X) {
	T x;
	HANDLE_ERROR(cudaMemcpy(&x,X,sizeof(T),cudaMemcpyDeviceToHost));
	return x;
}

template double* femib::cuda::copyToDevice<double>(double x);
template double femib::cuda::copyToHost<double>(double *x);

