find_package (benchmark REQUIRED)
find_package (Eigen3 REQUIRED)
find_package (TBB REQUIRED)

add_custom_target (bench COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS gauss_bench affine_bench cuda_bench)

include_directories (../src/gauss)
add_executable (gauss_bench gauss_bench.cpp)
target_compile_features (gauss_bench PRIVATE cxx_std_17)
target_link_libraries (gauss_bench Eigen3::Eigen TBB::tbb benchmark::benchmark)

include_directories (../src/affine)
add_executable (affine_bench affine_bench.cpp)
target_compile_features (affine_bench PRIVATE cxx_std_17)
target_link_libraries (affine_bench Eigen3::Eigen TBB::tbb benchmark::benchmark)

include_directories (../src/cuda)
add_executable (cuda_bench cuda_bench.cpp)
target_compile_definitions(cuda_bench PRIVATE MESH_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../mesh/")
target_compile_features (cuda_bench PRIVATE cxx_std_17)
target_link_libraries (cuda_bench cuda read Eigen3::Eigen TBB::tbb benchmark::benchmark)
