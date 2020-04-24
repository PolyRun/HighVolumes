# ---------------------- cmake setup

cmake_minimum_required(VERSION 3.6)

# --- add files to volume target in parent directory

add_subdirectory(dotProduct)

target_sources(volume PRIVATE
volume/linalg/cholesky.c
volume/linalg/linalg.c
volume/linalg/cholesky.h
volume/linalg/linalg.h)


