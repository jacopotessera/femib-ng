cmake_minimum_required(VERSION 3.9.1)
set(LIBRARY_OUTPUT_PATH  ${CMAKE_BINARY_DIR}/lib)
add_library(read SHARED read.cpp)
find_package (spdlog)
if(spdlog_FOUND)
       target_link_libraries (read spdlog::spdlog)
elseif(NOT spdlog_FOUND)
        message(FATAL_ERROR "
FATAL:
        spdlog Not Found.
	")
endif()
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
if(Eigen3_FOUND)
    target_link_libraries (read Eigen3::Eigen)
elseif(NOT Eigen3_FOUND)
        message(FATAL_ERROR "
FATAL:
        Eigen3 Not Found.
")
endif()
