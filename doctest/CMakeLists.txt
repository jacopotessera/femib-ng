# Prepare "Catch" library for other executables
find_package (doctest REQUIRED)
include_directories (${TEST_SOURCE_DIR}/src/affine)

add_executable (DocTest test.cpp)
target_compile_features(DocTest PRIVATE cxx_std_17)
target_link_libraries (DocTest
                       affine
		       doctest::doctest
		       )
	       add_test (NAME MyDocTest COMMAND DocTest)
	       #               add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND}
	       #       DEPENDS DocTest)
