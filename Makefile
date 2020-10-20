.PHONY: clean prepare cmake lib

CMAKE_LISTS_TXT := CMakeLists.txt
SRC_DIR := src
TEST_DIR := test
LIB_SRC_DIR := $(SRC_DIR)/$(LIB)
LIB_CMAKE := $(LIB_SRC_DIR)/$(CMAKE_LIST_TXT)
LIB_HPP := $(LIB_SRC_DIR)/$(LIB).hpp
LIB_CPP := $(LIB_SRC_DIR)/$(LIB).cpp
LIB_TEST_CPP := $(TEST_DIR)/$(LIB)_test.cpp
LIB_TEST_CMAKE := $(TEST_DIR)/$(CMAKE_LISTS_TXT)
LIB_UPPERCASE := $(shell echo $(LIB) | tr a-z A-Z)

clean:
	rm -rf build/

prepare:
	mkdir -p build

cmake: clean prepare
	cmake -H. -Bbuild

lib:
	#CMAKE
	echo 'add_subdirectory($${SRC_DIR}/$(LIB))' >> $(CMAKE_LISTS_TXT)

	mkdir -p $(LIB_SRC_DIR)
	#LIB CMAKE
	touch $(LIB_CMAKE)
	echo 'cmake_minimum_required(VERSION 3.9.1)' > src/$(LIB)/CMakeLists.txt
	echo 'set(LIBRARY_OUTPUT_PATH  $${CMAKE_BINARY_DIR}/lib)' >> src/$(LIB)/CMakeLists.txt
	echo 'add_library($(LIB) SHARED $(LIB).cpp)' >> src/$(LIB)/CMakeLists.txt
	echo 'find_package (spdlog)' >> src/$(LIB)/CMakeLists.txt
	echo 'if(spdlog_FOUND)' >> src/$(LIB)/CMakeLists.txt
	echo '   message("spdlog found.")' >> src/$(LIB)/CMakeLists.txt
	echo '       target_link_libraries ($(LIB) spdlog::spdlog)' >> src/$(LIB)/CMakeLists.txt
	echo 'elseif(NOT spdlog_FOUND)' >> src/$(LIB)/CMakeLists.txt
	echo '        message(FATAL_ERROR "' >> src/$(LIB)/CMakeLists.txt
	echo 'FATAL:' >> src/$(LIB)/CMakeLists.txt
	echo '        spdlog Not Found.' >> src/$(LIB)/CMakeLists.txt
	echo '	")' >> src/$(LIB)/CMakeLists.txt
	echo 'endif()' >> src/$(LIB)/CMakeLists.txt
	#LIB HPP
	touch $(LIB_HPP)
	echo '#ifndef $(LIB_UPPERCASE)_HPP_INCLUDED_' > $(LIB_HPP)
	echo '#define $(LIB_UPPERCASE)_HPP_INCLUDED_' >> $(LIB_HPP)
	echo 'namespace femib::$(LIB) {' >> $(LIB_HPP)
	echo '' >> $(LIB_HPP)
	echo '}' >> $(LIB_HPP)
	echo '#endif' >> $(LIB_HPP)
	#LIB CPP
	touch $(LIB_CPP)
	echo '#include "$(LIB).hpp"' > $(LIB_CPP)
	echo '#include "spdlog/spdlog.h"' > $(LIB_CPP)
	#TEST CMAKE
	sed -i '/add_custom_target /s/)$$//' $(LIB_TEST_CMAKE)
	sed -i '/add_custom_target /s/$$/ $(LIB)_test)/' $(LIB_TEST_CMAKE)
	echo '' >> $(LIB_TEST_CMAKE)
	echo 'include_directories (../src/$(LIB))' >> $(LIB_TEST_CMAKE)
	echo 'add_executable ($(LIB)_test $(LIB)_test.cpp)' >> $(LIB_TEST_CMAKE)
	echo 'target_compile_features ($(LIB)_test PRIVATE cxx_std_17)' >> $(LIB_TEST_CMAKE)
	echo 'target_link_libraries ($(LIB)_test $(LIB) doctest::doctest)' >> $(LIB_TEST_CMAKE)
	echo 'add_test (NAME $(LIB)_test COMMAND $(LIB)_test)' >> $(LIB_TEST_CMAKE)
	#TEST CPP
	touch $(LIB_TEST_CPP)
	echo '#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN' > $(LIB_TEST_CPP)
	echo '#include <doctest/doctest.h>' >> $(LIB_TEST_CPP)
	echo '#include "../$(LIB_HPP)"' >> $(LIB_TEST_CPP)
	echo '' >> $(LIB_TEST_CPP)
	echo 'TEST_CASE("testing $(LIB)") {' >> $(LIB_TEST_CPP)
	echo '    CHECK(true);' >> $(LIB_TEST_CPP)
	echo '}' >> $(LIB_TEST_CPP)
