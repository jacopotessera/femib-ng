cmake_minimum_required(VERSION 3.9.1)

enable_language(CUDA)
include(CheckLanguage)
check_language(CUDA)

if(NOT DEFINED CMAKE_CUDA_STANDARD)
    set(CMAKE_CUDA_STANDARD 11)
    set(CMAKE_CUDA_ARCHITECTURES 61)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
endif()

add_library(cuda SHARED cuda.cu)

find_package (spdlog)
if(spdlog_FOUND)
    target_link_libraries (cuda spdlog::spdlog)
elseif(NOT spdlog_FOUND)
	message(FATAL_ERROR "
FATAL:
	spdlog Not Found.
")
endif()

set(LIBRARY_OUTPUT_PATH  ${CMAKE_BINARY_DIR}/lib)

