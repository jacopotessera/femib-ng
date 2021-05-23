find_package (benchmark REQUIRED)

add_custom_target (bench COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS gauss_bench)

include_directories (../src/gauss)
add_executable (gauss_bench gauss_bench.cpp)
target_compile_features (gauss_bench PRIVATE cxx_std_17)
target_link_libraries (gauss_bench gauss benchmark::benchmark)
#add_test (NAME gauss_bench COMMAND gauss_bench)

