# ---------------------- cmake setup

cmake_minimum_required(VERSION 3.6)

# MB: maybe not needed?
# LINK_DIRECTORIES(/usr/include)


# ----------------------- prng
set(PRNG_SOURCE_FILES
   random/prng.c
)
set(PRNG_HEADER_FILES
   random/prng.h
)
add_library(prng ${PRNG_SOURCE_FILES} ${PRNG_HEADER_FILES})


# ----------------------- volume
set(VOLUME_SOURCE_FILES
   volume/volume.c volume/preprocess/init_ellipsoid.c volume/preprocess/preprocess.c volume/cholesky.c volume/matrix.c volume/polytope.c volume/ellipsoid.c volume/global.c
)
set(VOLUME_HEADER_FILES
   volume/volume.h volume/preprocess.h volume/cholesky.h volume/matrix.h volume/polytope.h volume/ellipsoid.h volume/global.h
)
add_library(volume ${VOLUME_SOURCE_FILES} ${VOLUME_HEADER_FILES})
target_link_libraries(volume PUBLIC glpk prng)

set_target_properties(volume PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
   COMPILE_FLAGS "-march=native -mfma -ffast-math -O3"
)

# ----------------------- volume_helper
set(VOLUME_HELPER_SOURCE_FILES
   volume/volume_helper.cpp
)
set(VOLUME_HELPER_HEADER_FILES
   volume/volume_helper.hpp
)
add_library(volume_helper ${VOLUME_HELPER_SOURCE_FILES} ${VOLUME_HELPER_HEADER_FILES})

target_link_libraries(volume_helper PUBLIC volume)

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

target_link_libraries(silaeubl PUBLIC prng)

set_target_properties(silaeubl PROPERTIES
   CXX_STANDARD 11
   CXX_STANDARD_REQUIRED YES
   CXX_EXTENSIONS NO
)



# ----------------------- test_preprocess
#project(test_preprocess)

#set(TEST_PREPROCESS_SOURCE_FILE
#   test_preprocess.cpp
#)
#add_executable(test_preprocess ${TEST_PREPROCESS_SOURCE_FILE})

#target_link_libraries(test_preprocess PUBLIC volume volume_helper)

#set_target_properties(test_preprocess PROPERTIES
#   CXX_STANDARD 11
#   CXX_STANDARD_REQUIRED YES
#   CXX_EXTENSIONS NO
#)
