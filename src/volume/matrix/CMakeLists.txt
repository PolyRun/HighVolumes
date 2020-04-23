# ---------------------- cmake setup

cmake_minimum_required(VERSION 3.6)

# ----------------------- matrix
set(MATRIX_SOURCE_FILES
   matrix.c
)
set(MATRIX_HEADER_FILES
   matrix.h
)
add_library(matrix ${MATRIX_SOURCE_FILES} ${MATRIX_HEADER_FILES})

