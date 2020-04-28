cmake_minimum_required(VERSION 3.10)


project(test_end_to_end)

add_executable(test_end_to_end end_to_end_test.cpp)

target_link_libraries(test_end_to_end PUBLIC volume_helper)

set_target_properties(test_end_to_end PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
)



