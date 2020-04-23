# ---------------------- cmake setup

cmake_minimum_required(VERSION 3.6)

# --- add files to volume target in parent directory

target_sources(volume PRIVATE
volume/polytope/polytope.c
volume/polytope/polytope_glpk.c
volume/polytope/polytope.h
)
