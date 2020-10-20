cmake_minimum_required(VERSION 3.9.1)
set(LIBRARY_OUTPUT_PATH  ${CMAKE_BINARY_DIR}/lib)
add_library(gauss SHARED gauss.cpp)
find_package (spdlog)
if(spdlog_FOUND)
   message("spdlog found.")
       target_link_libraries (gauss spdlog::spdlog)
elseif(NOT spdlog_FOUND)
        message(FATAL_ERROR "
FATAL:
        spdlog Not Found.
	")
endif()
