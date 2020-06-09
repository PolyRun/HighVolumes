cmake_minimum_required(VERSION 3.10)

file(COPY benchmark.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY bench_birk.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY plot.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY bm.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(build-time-make-directory ALL
   COMMAND ${CMAKE_COMMAND} -E make_directory out plots)

find_package(OpenMP)

project(benchmark_polyvest)

add_executable(benchmark_polyvest benchmark_polyvest.cpp)

target_link_libraries(benchmark_polyvest PUBLIC volume_helper polyvest_lib OpenMP::OpenMP_CXX)

set_target_properties(benchmark_polyvest PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
)

project(benchmark_A1_volume)

add_executable(benchmark_A1_volume benchmark_A1_volume.cpp)

target_link_libraries(benchmark_A1_volume PUBLIC volume_helper polyvest_lib OpenMP::OpenMP_CXX)

set_target_properties(benchmark_A1_volume PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
)

project(benchmark_A2_preprocessing)

add_executable(benchmark_A2_preprocessing benchmark_A2_preprocessing.cpp)

target_link_libraries(benchmark_A2_preprocessing PUBLIC volume_helper polyvest_lib OpenMP::OpenMP_CXX)

set_target_properties(benchmark_A2_preprocessing PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
)


project(benchmark_app)

add_executable(benchmark_app benchmark_app.cpp)

target_link_libraries(benchmark_app PUBLIC volume_helper polyvest_lib OpenMP::OpenMP_CXX)

set_target_properties(benchmark_app PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
)



project(benchmark_dotProduct)

add_executable(benchmark_dotProduct benchmark_dotProduct.cpp)

target_link_libraries(benchmark_dotProduct PUBLIC volume_helper OpenMP::OpenMP_CXX)

set_target_properties(benchmark_dotProduct PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
)


project(benchmark_intersect)

add_executable(benchmark_intersect benchmark_intersect.cpp)

target_link_libraries(benchmark_intersect PUBLIC volume_helper OpenMP::OpenMP_CXX)

set_target_properties(benchmark_intersect PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
   COMPILE_FLAGS "-march=native -mfma -ffast-math -O3"
)


project(benchmark_jit)

add_executable(benchmark_jit benchmark_jit.cpp)

target_link_libraries(benchmark_jit PUBLIC volume_helper OpenMP::OpenMP_CXX)

set_target_properties(benchmark_jit PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
)



project(benchmark_ball_intersect)

add_executable(benchmark_ball_intersect benchmark_ball_intersect.cpp)

target_link_libraries(benchmark_ball_intersect PUBLIC volume_helper OpenMP::OpenMP_CXX)

set_target_properties(benchmark_ball_intersect PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
)



project(benchmark_preprocess)

add_executable(benchmark_preprocess benchmark_preprocessing.cpp)

target_link_libraries(benchmark_preprocess PUBLIC volume_helper OpenMP::OpenMP_CXX)

set_target_properties(benchmark_preprocess PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
)


project(benchmark_A1_polytope)

add_executable(benchmark_A1_polytope benchmark_A1_polytope.cpp)

target_link_libraries(benchmark_A1_polytope PUBLIC volume_helper OpenMP::OpenMP_CXX)

set_target_properties(benchmark_A1_polytope PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
)


project(convergence)

add_executable(convergence convergence.cpp)

target_link_libraries(convergence PUBLIC volume_helper OpenMP::OpenMP_CXX)

set_target_properties(convergence PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
)

project(benchmark_randomness)

add_executable(benchmark_randomness benchmark_randomness.cpp)

target_link_libraries(benchmark_randomness PUBLIC volume_helper prng OpenMP::OpenMP_CXX)

set_target_properties(benchmark_randomness PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
)
