#include <stdio.h>

static void HandleError(cudaError_t err, const char *file, int line )
{
    if(err != cudaSuccess)
        {
        printf( "%s in %s at line %d\n", cudaGetErrorString( err ), file, line );
        exit( EXIT_FAILURE );
    }
}

#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))


#include <iostream>

#include "cuda.h"

void femib::cuda::printSize() {
	std::cout << "[CUDA stack size] found to be " << getStackSize() << std::endl;
	std::cout << "[CUDA heap  size] found to be " << getHeapSize()  << std::endl;
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

