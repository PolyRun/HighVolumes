# LOG for VTUNE images

## 1 CSC Sparse
benchmark/benchmark_A1_volume -b "generator=2var_100,polytopeOptimize=true,polytopeType=2,r=1" -c "step_size=400" -f "PolytopeCSC_intersectCoord=cached_b_vec_nan_inv,rand_f=sr_rand"

## 2 CSC Sparse
benchmark/benchmark_A1_volume -b "generator=2var_100,polytopeOptimize=true,polytopeType=2,r=1" -c "step_size=400" -f "PolytopeCSC_intersectCoord=cached_b_vec_nan_inv,rand_f=std_rand"

## 3 JIT Sparse
benchmark/benchmark_A1_volume -b "generator=2var_100,polytopeOptimize=true,polytopeType=3,r=1" -c "step_size=400" -f "PolytopeJIT_gen=quad_data_acc,rand_f=sr_rand"

## 4 JIT Sparse
benchmark/benchmark_A1_volume -b "generator=2var_100,polytopeOptimize=true,polytopeType=3,r=1" -c "step_size=400" -f "PolytopeJIT_gen=quad_data_acc,rand_f=std_rand"

Note: it is clear: we have fixed the random issue. Now other issues are more important!


## 5 JIT Dense
-b "generator=cube_rot_r1.0_60,polytopeOptimize=true,polytopeType=3,r=1" -c "step_size=1000" -f "PolytopeJIT_gen=quad_data_acc,rand_f=sr_rand"

I assume the 23% "unknown stack frames" are mostly due to the JIT code that cannot be linked to any source files.

## 6 CSC Dense
-b "generator=cube_rot_r1.0_100,polytopeOptimize=true,polytopeType=2,r=1" -c "step_size=1000" -f "PolytopeCSC_intersectCoord=cached_b_vec_nan_inv,rand_f=sr_rand"

## 7 PolytopeT
-b "generator=cube_rot_r1.0_100,polytopeOptimize=true,polytopeType=0,r=1" -c "step_size=1000" -f "PolytopeT_intersectCoord=cached_b_vec,rand_f=sr_rand"

## 8 Ellipsoid
TODO


