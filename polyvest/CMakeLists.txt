cmake_minimum_required(VERSION 3.6)

project(polyvest_lib)

set(POLYVESTLIB_SOURCE_FILES
   vol.cpp
)
set(POLYVESTLIB_HEADER_FILES
   vol.h
)

add_library(polyvest_lib ${POLYVESTLIB_SOURCE_FILES} ${POLYVESTLIB_HEADER_FILES})
target_link_libraries(polyvest_lib PUBLIC glpk armadillo)



project(polyvest)

add_executable(polyvest main.cpp)

target_link_libraries(polyvest PUBLIC polyvest_lib)

set_target_properties(polyvest PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
)
