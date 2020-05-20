cmake_minimum_required(VERSION 3.10)

add_library(test_helper test_helpers.cpp test_helpers.hpp)

add_subdirectory(general)
add_subdirectory(volume)
add_subdirectory(preprocess)
add_subdirectory(matrix)
add_subdirectory(polytopeCSC)
add_subdirectory(jit)
add_subdirectory(random)


file(COPY test.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
