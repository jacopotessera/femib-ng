cmake_minimum_required(VERSION 3.9.1)
set(LIBRARY_OUTPUT_PATH  ${CMAKE_BINARY_DIR}/lib)
add_library(mongo SHARED mongo.cpp)
find_package (spdlog)
if(spdlog_FOUND)
       target_link_libraries (mongo spdlog::spdlog)
elseif(NOT spdlog_FOUND)
        message(FATAL_ERROR "
FATAL:
        spdlog Not Found.
")
endif()
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
if(Eigen3_FOUND)
target_link_libraries (mongo Eigen3::Eigen)
elseif(NOT Eigen3_FOUND)
        message(FATAL_ERROR "
FATAL:
        Eigen3 Not Found.
")
endif()

find_package(libmongocxx REQUIRED)
find_package(libbsoncxx REQUIRED)
include_directories(${LIBMONGOCXX_INCLUDE_DIR})
include_directories(${LIBBSONCXX_INCLUDE_DIR})
include_directories("/usr/include/mongocxx/v_noabi")
include_directories("/usr/include/bsoncxx/v_noabi")

target_link_libraries(mongo ${LIBMONGOCXX_LIBRARIES})
target_link_libraries(mongo ${LIBBSONCXX_LIBRARIES})

