cmake_minimum_required(VERSION 3.10)

file(COPY benchmark.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY plot.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(build-time-make-directory ALL
   COMMAND ${CMAKE_COMMAND} -E make_directory out plots)

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

project(benchmark_A1_volume)

add_executable(benchmark_A1_volume benchmark_A1_volume.cpp)

target_link_libraries(benchmark_A1_volume PUBLIC volume_helper)

set_target_properties(benchmark_A1_volume PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
)


project(benchmark_dotProduct)

add_executable(benchmark_dotProduct benchmark_dotProduct.cpp)

target_link_libraries(benchmark_dotProduct PUBLIC volume_helper)

set_target_properties(benchmark_dotProduct PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
)


project(benchmark_intersect)

add_executable(benchmark_intersect benchmark_intersect.cpp)

target_link_libraries(benchmark_intersect PUBLIC volume_helper)

set_target_properties(benchmark_intersect PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
)


project(benchmark_preprocess)

add_executable(benchmark_preprocess benchmark_preprocessing.cpp)

target_link_libraries(benchmark_preprocess PUBLIC volume_helper)

set_target_properties(benchmark_preprocess PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
)


project(benchmark_A1_polytope)

add_executable(benchmark_A1_polytope benchmark_A1_polytope.cpp)

target_link_libraries(benchmark_A1_polytope PUBLIC volume_helper)

set_target_properties(benchmark_A1_polytope PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
)
