find_package (benchmark REQUIRED)
find_package (Eigen3 REQUIRED)

add_custom_target (bench COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS gauss_bench)

include_directories (../src/gauss)
add_executable (gauss_bench gauss_bench.cpp)
target_compile_features (gauss_bench PRIVATE cxx_std_17)
target_link_libraries (gauss_bench gauss benchmark::benchmark)

include_directories (../src/affine)
add_executable (affine_bench affine_bench.cpp)
target_compile_features (affine_bench PRIVATE cxx_std_17)
target_link_libraries (affine_bench affine Eigen3::Eigen benchmark::benchmark)
