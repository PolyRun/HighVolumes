# ---------------------- cmake setup

cmake_minimum_required(VERSION 3.6)

# MB: maybe not needed?
# LINK_DIRECTORIES(/usr/include)

# ----------------------- volume
set(VOLUME_SOURCE_FILES
   poly/volume.c beta_cut.c poly/cube.c
)
set(VOLUME_HEADER_FILES
   poly/volume.h beta_cut.h poly/cube.h
)
add_library(volume ${VOLUME_SOURCE_FILES} ${VOLUME_HEADER_FILES})

# ----------------------- peterem
project(peterem)

set(PETEREM_HEADER_FILES
   peterem.hpp
)
set(PETEREM_SOURCE_FILE
   peterem.cpp
)
add_executable(peterem ${PETEREM_SOURCE_FILE} ${PETEREM_HEADER_FILES})

target_link_libraries(peterem PUBLIC volume glpk)

set_target_properties(peterem PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
)

