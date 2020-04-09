cmake_minimum_required(VERSION 3.10)

add_subdirectory(volume)
add_subdirectory(preprocess)

file(COPY test.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
