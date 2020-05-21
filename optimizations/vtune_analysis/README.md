# LOG for VTUNE images

## 1
benchmark/benchmark_A1_volume -b "generator=2var_100,polytopeOptimize=true,polytopeType=2,r=1" -c "step_size=400" -f "PolytopeCSC_intersectCoord=cached_b_vec_nan_inv,rand_f=sr_rand"

## 2
benchmark/benchmark_A1_volume -b "generator=2var_100,polytopeOptimize=true,polytopeType=2,r=1" -c "step_size=400" -f "PolytopeCSC_intersectCoord=cached_b_vec_nan_inv,rand_f=std_rand"

## 3
benchmark/benchmark_A1_volume -b "generator=2var_100,polytopeOptimize=true,polytopeType=3,r=1" -c "step_size=400" -f "PolytopeJIT_gen=quad_data_acc,rand_f=sr_rand"

## 4
benchmark/benchmark_A1_volume -b "generator=2var_100,polytopeOptimize=true,polytopeType=3,r=1" -c "step_size=400" -f "PolytopeJIT_gen=quad_data_acc,rand_f=std_rand"



