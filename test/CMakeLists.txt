find_package (Boost COMPONENTS system filesystem unit_test_framework REQUIRED)
set(TEST_SOURCE_DIR ..)
message(${TEST_SOURCE_DIR}/src/affine)
include_directories (${TEST_SOURCE_DIR}/src/affine
                     ${Boost_INCLUDE_DIRS}
                     )
add_definitions (-DBOOST_TEST_DYN_LINK)
add_executable (Test test.cpp)
target_link_libraries (Test
                       affine
                       ${Boost_FILESYSTEM_LIBRARY}
                       ${Boost_SYSTEM_LIBRARY}
                       ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
                       )
	       add_test (NAME MyTest COMMAND Test)
	       add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND}
                  DEPENDS Test)
