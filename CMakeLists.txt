cmake_minimum_required(VERSION 3.9.1)

project(femib LANGUAGES CXX CUDA)

set(SRC_DIR src)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Disable in-source builds to prevent source tree corruption.
if(" ${CMAKE_SOURCE_DIR}" STREQUAL " ${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "
FATAL: In-source builds are not allowed.
       You should create a separate directory for build files.
")
endif()

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH  ${CMAKE_BINARY_DIR}/lib)

SET(GCC_COVERAGE_COMPILE_FLAGS "-fprofile-arcs -ftest-coverage --coverage -g")
SET(GCC_COVERAGE_LINK_FLAGS    "-lgcov")
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")
SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}")

add_executable(femib ${SRC_DIR}/main.cpp)
target_link_libraries(femib)

include(CodeCoverage)
append_coverage_compiler_flags()
setup_target_for_coverage_gcovr_html(NAME coverage COMMAND check DEPENDENCIES check BASE_DIRECTORY src/)

enable_testing ()
add_subdirectory(test)

add_subdirectory(${SRC_DIR}/affine)
add_subdirectory(${SRC_DIR}/cuda)
add_subdirectory(${SRC_DIR}/gauss)
add_subdirectory(${SRC_DIR}/types)
add_subdirectory(${SRC_DIR}/mesh)
add_subdirectory(${SRC_DIR}/finite_element)
add_subdirectory(${SRC_DIR}/finite_element_space)
