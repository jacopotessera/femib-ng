.PHONY: clean prepare cmake lib

clean:
	rm -rf build/

prepare:
	mkdir -p build

cmake: clean prepare
	cmake -H. -Bbuild

lib:
	mkdir -p src/$(LIB)
	touch src/$(LIB)/CMakeLists.txt
	echo 'cmake_minimum_required(VERSION 3.9.1)' > src/$(LIB)/CMakeLists.txt
	echo 'set(LIBRARY_OUTPUT_PATH  ${CMAKE_BINARY_DIR}/lib)' >> src/$(LIB)/CMakeLists.txt
	echo 'add_library($(LIB) SHARED $(LIB).cpp)' >> src/$(LIB)/CMakeLists.txt
	echo 'find_package (spdlog)' >> src/$(LIB)/CMakeLists.txt
	echo 'if(spdlog_FOUND)' >> src/$(LIB)/CMakeLists.txt
	echo '   message("spdlog found.")' >> src/$(LIB)/CMakeLists.txt
	echo '       target_link_libraries (cuda spdlog::spdlog)' >> src/$(LIB)/CMakeLists.txt
	echo 'elseif(NOT spdlog_FOUND)' >> src/$(LIB)/CMakeLists.txt
	echo '        message(FATAL_ERROR "' >> src/$(LIB)/CMakeLists.txt
	echo 'FATAL:' >> src/$(LIB)/CMakeLists.txt
	echo '        spdlog Not Found.' >> src/$(LIB)/CMakeLists.txt
	echo '	")' >> src/$(LIB)/CMakeLists.txt
	echo 'endif()' >> src/$(LIB)/CMakeLists.txt
	echo '' >> src/$(LIB)/CMakeLists.txt
	touch src/$(LIB)/$(LIB).cpp
	touch src/$(LIB)/$(LIB).hpp
