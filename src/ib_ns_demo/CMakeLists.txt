find_package (spdlog REQUIRED)
find_package (TBB REQUIRED)
find_package (Eigen3 REQUIRED NO_MODULE)

add_executable (ib_ns_demo ib_ns_demo.cpp)
target_compile_definitions(ib_ns_demo PRIVATE MESH_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../../mesh/")
target_compile_features (ib_ns_demo PRIVATE cxx_std_23)
target_link_libraries (ib_ns_demo read cuda spdlog::spdlog Eigen3::Eigen TBB::tbb)
