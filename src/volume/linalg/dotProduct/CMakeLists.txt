# ---------------------- cmake setup

cmake_minimum_required(VERSION 3.6)

# --- add files to volume target in parent directory

target_sources(volume PRIVATE
volume/linalg/dotProduct/dotProduct.h
volume/linalg/dotProduct/dotProduct.c
volume/linalg/dotProduct/dotProduct_auto1.c
volume/linalg/dotProduct/dotProduct_auto2.c
)


