cmake_minimum_required(VERSION 2.8)
set(CMAKE_CXX_FLAGS "-std=c++11 -Werror -Wall -Wfatal-errors")

add_executable(test_simple test_simple.cpp kdarray.hpp)
