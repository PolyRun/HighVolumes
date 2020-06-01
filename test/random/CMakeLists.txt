cmake_minimum_required(VERSION 3.10)

project(test_random)

add_executable(test_random random.cpp)

target_link_libraries(test_random PUBLIC prng volume_helper)

set_target_properties(test_random PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
   COMPILE_FLAGS "-march=native"
)


