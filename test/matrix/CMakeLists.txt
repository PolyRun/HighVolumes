cmake_minimum_required(VERSION 3.10)

project(matrix_unit_test)

add_executable(matrix_unit_test unittests.cpp)

target_link_libraries(matrix_unit_test PUBLIC volume_helper test_helper)


set_target_properties(matrix_unit_test PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
)
