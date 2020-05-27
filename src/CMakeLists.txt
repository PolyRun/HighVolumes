# ---------------------- cmake setup

cmake_minimum_required(VERSION 3.6)

find_package(OpenMP)

# ----------------------- prng
set(PRNG_SOURCE_FILES
   random/prng.c
)
set(PRNG_HEADER_FILES
   random/prng.h
)
add_library(prng ${PRNG_SOURCE_FILES} ${PRNG_HEADER_FILES})

set_target_properties(prng PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
   COMPILE_FLAGS "-march=native -mavx -mfma -ffast-math -O3"
)


# ----------------------- jit
set(JIT_SOURCE_FILES
   jit/jit.c
)
set(JIT_HEADER_FILES
   jit/jit.h
)
add_library(jit ${JIT_SOURCE_FILES} ${JIT_HEADER_FILES})


# ----------------------- volume



set(VOLUME_SOURCE_FILES
   volume/volume.c volume/preprocess/init_ellipsoid.c volume/preprocess/preprocess.c 
)
set(VOLUME_HEADER_FILES
   volume/volume.h volume/preprocess/preprocess.h
)

add_library(volume ${VOLUME_SOURCE_FILES} ${VOLUME_HEADER_FILES})

add_subdirectory(volume/linalg)
add_subdirectory(volume/matrix)
add_subdirectory(volume/ellipsoid)
add_subdirectory(volume/polytope)
add_subdirectory(volume/polytopeT)
add_subdirectory(volume/polytopeCSC)
add_subdirectory(volume/polytopeJIT)

target_link_libraries(volume PRIVATE glpk prng jit OpenMP::OpenMP_CXX)

set_target_properties(volume PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
   COMPILE_FLAGS "-march=native -mfma -ffast-math -O3"
)

# ----------------------- performance_counter
set(PERFORMANCE_COUNTER_SOURCES_FILES
   util/performance_counter.cpp
)
set(PERFORMANCE_COUNTER_HEADER_FILES
   util/performance_counter.hpp
)
add_library(performance_counter ${PERFORMANCE_COUNTER_SOURCES_FILES} ${PERFORMANCE_COUNTER_HEADER_FILES})

set_target_properties(performance_counter PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
)

# ----------------------- volume_helper
set(VOLUME_HELPER_SOURCE_FILES
   volume/volume_helper.cpp volume/volume_cost.cpp volume/volume_examples.cpp
)
set(VOLUME_HELPER_HEADER_FILES
   volume/volume_helper.hpp
)
add_library(volume_helper ${VOLUME_HELPER_SOURCE_FILES} ${VOLUME_HELPER_HEADER_FILES})

target_link_libraries(volume_helper PUBLIC volume performance_counter)

set_target_properties(volume_helper PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
)

# ----------------------- peterem
project(peterem)

set(PETEREM_HEADER_FILES
   peterem.hpp
)
set(PETEREM_SOURCE_FILE
   peterem.cpp
)
add_executable(peterem ${PETEREM_SOURCE_FILE} ${PETEREM_HEADER_FILES})

target_link_libraries(peterem PUBLIC volume_helper)

set_target_properties(peterem PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
)


# ----------------------- silaeubl
project(silaeubl)

set(SILAEUBL_HEADER_FILES
   silaeubl.hpp
)
set(SILAEUBL_SOURCE_FILE
   silaeubl.cpp
)
add_executable(silaeubl ${SILAEUBL_SOURCE_FILE} ${SILAEUBL_HEADER_FILES})

target_link_libraries(silaeubl PUBLIC volume_helper prng)

add_subdirectory(random)

set_target_properties(silaeubl PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
)


# ---------------------- random poly to file
add_subdirectory(util)
