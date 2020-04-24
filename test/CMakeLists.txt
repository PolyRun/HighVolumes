cmake_minimum_required(VERSION 3.10)

add_library(test_helper test_helpers.cpp test_helpers.hpp)

add_subdirectory(volume)
add_subdirectory(preprocess)
add_subdirectory(matrix)


file(COPY test.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
