# ---------------------- cmake setup

cmake_minimum_required(VERSION 3.6)

# ----------------------- matrix
set(ELLIPSOID_SOURCE_FILES
   ellipsoid.c
)
set(ELLIPSOID_HEADER_FILES
   ellipsoid.h
)
add_library(ellipsoid ${ELLIPSOID_SOURCE_FILES} ${ELLIPSOID_HEADER_FILES})

