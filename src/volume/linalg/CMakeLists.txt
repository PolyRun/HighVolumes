# ---------------------- cmake setup

cmake_minimum_required(VERSION 3.6)

# --- add files to volume target in parent directory

target_sources(volume PUBLIC
${CMAKE_CURRENT_SOURCE_DIR}/cholesky.c
${CMAKE_CURRENT_SOURCE_DIR}/linalg.c
${CMAKE_CURRENT_SOURCE_DIR}/cholesky.h
${CMAKE_CURRENT_SOURCE_DIR}/linalg.h)


