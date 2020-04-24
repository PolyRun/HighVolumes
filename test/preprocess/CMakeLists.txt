cmake_minimum_required(VERSION 3.10)

project(test_preprocess)

OPTION(OUTPUT_MACRO "all debug output in source files is printed" OFF)

IF(OUTPUT_MACRO)
    ADD_DEFINITIONS(-DDEBUG_MSG)
ENDIF(OUTPUT_MACRO)

add_executable(test_preprocess test_preprocess.cpp)

target_link_libraries(test_preprocess PUBLIC test_helper volume_helper polyvest_lib)

set_target_properties(test_preprocess PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
)



project(test_init)

add_executable(test_init test_init.cpp)

target_link_libraries(test_init PUBLIC test_helper volume_helper polyvest_lib)

set_target_properties(test_preprocess PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
)


project(unit_test)

add_executable(unit_test unittests.cpp)

target_link_libraries(unit_test PUBLIC test_helper volume_helper)


set_target_properties(unit_test PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
)
