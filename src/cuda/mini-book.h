#ifndef MINI_BOOK_H_INCLUDED_
#define MINI_BOOK_H_INCLUDED_

#include <spdlog/spdlog.h>

static void HandleError(cudaError_t err, const char *file, int line) {
  if (err != cudaSuccess) {
    spdlog::error("{} in {} at line {}", cudaGetErrorString(err), file, line);
    exit(EXIT_FAILURE);
  }
}

#define HANDLE_ERROR(err) (HandleError(err, __FILE__, __LINE__))

#endif
