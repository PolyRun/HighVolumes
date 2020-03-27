cmake_minimum_required(VERSION 3.6)


project(testvol)

add_executable(testvol test_preprocess.cpp)

target_link_libraries(testvol PUBLIC volume_helper polyvest_lib)

set_target_properties(testvol PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
)
