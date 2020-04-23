# ---------------------- cmake setup

cmake_minimum_required(VERSION 3.6)

# ----------------------- matrix
set(POLYTOPE_SOURCE_FILES
   polytopeT.c
)
set(ELLIPSOID_HEADER_FILES
   polytopeT.h
)
add_library(polytopeT ${POLYTOPE_SOURCE_FILES} ${POLYTOPE_HEADER_FILES})

