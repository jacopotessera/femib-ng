cmake_minimum_required(VERSION 3.9.1)

set(LIBRARY_OUTPUT_PATH  ${CMAKE_BINARY_DIR}/lib)

add_library(affine SHARED affine.cpp)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
# Check for libray, if found print message, include dirs and link libraries.
if(Eigen3_FOUND)
    message("Eigen3 found.")
    target_link_libraries (affine Eigen3::Eigen)
elseif(NOT Eigen3_FOUND)
	message(FATAL_ERROR "
FATAL:
	Eigen3 Not Found.
")
endif()

