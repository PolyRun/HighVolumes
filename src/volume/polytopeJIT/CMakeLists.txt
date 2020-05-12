# ---------------------- cmake setup

cmake_minimum_required(VERSION 3.6)

# --- add files to volume target in parent directory

target_sources(volume PRIVATE
	volume/polytopeJIT/polytopeJIT.c
	volume/polytopeJIT/polytopeJIT.h
)

add_subdirectory(jit)
