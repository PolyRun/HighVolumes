cmake_minimum_required(VERSION 3.10)

project(test_volume_basics)

add_executable(test_volume_basics basics.cpp)

target_link_libraries(test_volume_basics PUBLIC volume_helper)

set_target_properties(test_volume_basics PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
   COMPILE_FLAGS "-march=native"
)

project(test_volume_estimate)

add_executable(test_volume_estimate estimate.cpp)

target_link_libraries(test_volume_estimate PUBLIC volume_helper)

set_target_properties(test_volume_estimate PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
)


project(test_volume_examples)

add_executable(test_volume_examples examples.cpp)

target_link_libraries(test_volume_examples PUBLIC volume_helper)

set_target_properties(test_volume_examples PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
)




