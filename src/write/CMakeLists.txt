cmake_minimum_required(VERSION 3.10)
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
add_library(write SHARED write.cpp)
target_link_libraries(write HighFive)
find_package (Eigen3 REQUIRED NO_MODULE)
if(Eigen3_FOUND)
    target_link_libraries (write Eigen3::Eigen)
elseif(NOT Eigen3_FOUND)
    message(FATAL_ERROR "
FATAL:
        Eigen3 Not Found.
")
endif()