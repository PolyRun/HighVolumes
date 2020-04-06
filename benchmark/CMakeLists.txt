cmake_minimum_required(VERSION 3.10)

project(benchmark)

add_executable(benchmark benchmark.cpp)

target_link_libraries(benchmark PUBLIC polyvest_lib volume_helper prng)

set_target_properties(benchmark PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
)



