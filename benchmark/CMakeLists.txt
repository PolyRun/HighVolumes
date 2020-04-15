cmake_minimum_required(VERSION 3.10)

add_custom_target(build-time-make-directory ALL
   COMMAND ${CMAKE_COMMAND} -E make_directory out)

project(benchmark_test_macro)

add_executable(benchmark_test_macro benchmark_test_macro.cpp)

target_link_libraries(benchmark_test_macro PUBLIC volume_helper)

set_target_properties(benchmark_test_macro PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
)



project(benchmark_test_xyz_f)

add_executable(benchmark_test_xyz_f benchmark_test_xyz_f.cpp)

target_link_libraries(benchmark_test_xyz_f PUBLIC volume_helper)

set_target_properties(benchmark_test_xyz_f PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
)



project(benchmark_polyvest)

add_executable(benchmark_polyvest benchmark_polyvest.cpp)

target_link_libraries(benchmark_polyvest PUBLIC volume_helper polyvest_lib)

set_target_properties(benchmark_polyvest PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
)