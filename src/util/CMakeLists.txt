cmake_minimum_required(VERSION 3.6)

project(random2file)

add_executable(random2file poly2file.cpp)

target_link_libraries(random2file PUBLIC volume_helper)

set_target_properties(peterem PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
)
