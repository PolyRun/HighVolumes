cmake_minimum_required(VERSION 3.10)

project(polytopecsc_unit_test)

add_executable(polytopecsc_unit_test unittests.cpp)

target_link_libraries(polytopecsc_unit_test PUBLIC volume_helper)


set_target_properties(polytopecsc_unit_test PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
)
