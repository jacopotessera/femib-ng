cmake_minimum_required(VERSION 3.10)
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
add_library(write SHARED write.cpp)
target_link_libraries(write HighFive)
