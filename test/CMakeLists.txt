# Prepare "Catch" library for other executables
find_package (doctest REQUIRED)

add_custom_target (check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS affine_test cuda_test)

include_directories (../src/affine)
add_executable (affine_test affine_test.cpp)
target_compile_features (affine_test PRIVATE cxx_std_17)
target_link_libraries (affine_test affine doctest::doctest)
add_test (NAME affine_test COMMAND affine_test)

include_directories (../src/cuda)
add_executable (cuda_test cuda_test.cpp)
target_compile_features (cuda_test PRIVATE cxx_std_17)
target_link_libraries (cuda_test cuda doctest::doctest)
add_test (NAME cuda_test COMMAND cuda_test)

