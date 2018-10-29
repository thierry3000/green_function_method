project(green_function_method)

cmake_minimum_required(VERSION 3.10)

if(CMAKE_BUILD_TYPE MATCHES Debug)
	message("debug output activated")
	set(CMAKE_CXX_FLAGS "-DDEBUG")
endif(CMAKE_BUILD_TYPE MATCHES Debug)

add_subdirectory(src/)
