# ---------------------- cmake setup

cmake_minimum_required(VERSION 3.6)

# --- add files to prng target in parent directory

target_sources(prng PRIVATE
random/generators/std_rand.h
random/generators/std_rand.c
random/generators/mersenne.h
random/generators/mersenne.c
random/generators/shiftreg.h
random/generators/shiftreg.c
)

