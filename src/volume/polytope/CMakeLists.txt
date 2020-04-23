# ---------------------- cmake setup

cmake_minimum_required(VERSION 3.6)

# ----------------------- matrix
set(POLYTOPE_SOURCE_FILES
   polytope.c
)
set(ELLIPSOID_HEADER_FILES
   polytope.h
)
add_library(polytope ${POLYTOPE_SOURCE_FILES} ${POLYTOPE_HEADER_FILES})

