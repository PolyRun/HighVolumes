cmake_minimum_required(VERSION 3.6)


project(test)

add_executable(test test_preprocess.cpp)

target_link_libraries(test PUBLIC volume_helper polyvest_lib)

set_target_properties(test PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
)
