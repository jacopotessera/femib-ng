cmake_minimum_required(VERSION 3.9.1)

project(femib LANGUAGES CUDA CXX)

set(SRC_DIR src)

set(CMAKE_CXX_STANDARD 14)

# Disable in-source builds to prevent source tree corruption.
if(" ${CMAKE_SOURCE_DIR}" STREQUAL " ${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "
FATAL: In-source builds are not allowed.
       You should create a separate directory for build files.
")
endif()

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH  ${CMAKE_BINARY_DIR}/lib)


add_executable(femib ${SRC_DIR}/main.cpp)
target_link_libraries(femib affine cuda)

enable_testing ()
add_subdirectory(test)

add_subdirectory(${SRC_DIR}/affine)
add_subdirectory(${SRC_DIR}/cuda)
add_subdirectory(${SRC_DIR}/gauss)
