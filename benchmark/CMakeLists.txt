cmake_minimum_required(VERSION 3.10)

project(test_macro)

add_executable(test_macro test_macro.cpp)

target_link_libraries(test_macro PUBLIC volume_helper)

set_target_properties(test_macro PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
)



project(test_mini)

add_executable(test_xyz_f test_xyz_f.cpp)

target_link_libraries(test_xyz_f PUBLIC volume_helper)

set_target_properties(test_xyz_f PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
)



