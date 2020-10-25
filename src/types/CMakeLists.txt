cmake_minimum_required(VERSION 3.9.1)
set(LIBRARY_OUTPUT_PATH  ${CMAKE_BINARY_DIR}/lib)
add_library(types SHARED types.cpp)
find_package (spdlog)
if(spdlog_FOUND)
   message("spdlog found.")
       target_link_libraries (types spdlog::spdlog)
elseif(NOT spdlog_FOUND)
        message(FATAL_ERROR "
FATAL:
        spdlog Not Found.
	")
endif()
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
# Check for libray, if found print message, include dirs and link libraries.
if(Eigen3_FOUND)
    message("Eigen3 found.")
    target_link_libraries (types Eigen3::Eigen)
elseif(NOT Eigen3_FOUND)
        message(FATAL_ERROR "
FATAL:
        Eigen3 Not Found.
")
endif()
