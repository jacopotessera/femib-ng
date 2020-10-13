#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include "cuda.h"
#include <iostream>

TEST_CASE("testing cuda") {
	femib::cuda::setStackSize(FEMIB_CUDA_STACK_SIZE);
	femib::cuda::setHeapSize(FEMIB_CUDA_HEAP_SIZE);
	femib::cuda::printSize();
	WARN(femib::cuda::getStackSize() == FEMIB_CUDA_STACK_SIZE);
	WARN(femib::cuda::getHeapSize() == FEMIB_CUDA_HEAP_SIZE);
	
	double x = 10;
	double *X = femib::cuda::copyToDevice<double>(x);
	double y = femib::cuda::copyToHost<double>(X);
	CHECK(x == y);
}

