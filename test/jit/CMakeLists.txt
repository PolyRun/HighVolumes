cmake_minimum_required(VERSION 3.10)

project(test_jit)

add_executable(test_jit jit.cpp)

target_link_libraries(test_jit PUBLIC jit)

set_target_properties(test_jit PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
   COMPILE_FLAGS "-march=native"
)


