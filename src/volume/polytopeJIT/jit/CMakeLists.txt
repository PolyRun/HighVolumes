# ---------------------- cmake setup

cmake_minimum_required(VERSION 3.6)

# --- add files to volume target in parent directory

target_sources(volume PRIVATE
	volume/polytopeJIT/jit/inside.h
	volume/polytopeJIT/jit/inside.c
	volume/polytopeJIT/jit/intersect.h
	volume/polytopeJIT/jit/intersect.c
	volume/polytopeJIT/jit/cacheReset.h
	volume/polytopeJIT/jit/cacheReset.c
	volume/polytopeJIT/jit/cacheUpdateCoord.h
	volume/polytopeJIT/jit/cacheUpdateCoord.c
	volume/polytopeJIT/jit/intersectCoord.h
	volume/polytopeJIT/jit/intersectCoord.c
)

