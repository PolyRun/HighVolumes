# ---------------------- cmake setup

cmake_minimum_required(VERSION 3.6)

# ----------------------- linalg
set(LINALG_SOURCE_FILES
   cholesky.c linalg.c
)
set(LINALG_HEADER_FILES
   cholesky.h linalg.h
)
add_library(linalg ${LINALG_SOURCE_FILES} ${LINALG_HEADER_FILES})

