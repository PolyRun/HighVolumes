
# turn off turboboost
# set memorybound to your streambenchmark value

export ERRBAR='On'

# 1) table
./benchmark/benchmark.py polytopeT_intersectCoord > polytopeT_intersectCoord.out

# 2) turn on errorbars
## python benchmark/benchmark.py polytopeT_cubeRot_bm
## python benchmark/benchmark.py polytopeT_cross_bm

# 3) maybe increase step_size from 100 to 1000 (line 726)
## ./benchmark/benchmark.py allbest_roofline_cube_rot_bm
## ./benchmark/benchmark.py allbest_roofline_cross_bm

# 4)
## python benchmark/benchmark.py ellipsoid_cacheUpdateCoord

# 6) turn on errorbars
## python benchmark/benchmark.py randomness

# 8) turn on errorbars
## python benchmark/benchmark.py csc_jit_only_cacheb_2var_intersectCoord_only


# 5) turn off errorbars! why? inherently large errors because of random element distribution in sparse matrix.
# TODO: add for n=200
export ERRBAR='Off'

export XLOG='On'
## python benchmark/benchmark.py density_runtime_bm
## python benchmark/benchmark.py density8_runtime_bm
export XLOG='Off'

# 7) turn off errorbars, maybe increase step_size from 100 to 1000 (line 659) or do multiple runs (line 661)
## ./benchmark/benchmark.py 2var_runtime_bm > 2var_runtime_bm.out
