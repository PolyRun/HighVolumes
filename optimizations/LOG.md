# Experiments Log

## Cross

```
benchmark/benchmark_A1_volume -b "r=1,polytopeTranspose=false,generator=cross_rn_8" -c "step_size=1000000"
``` 
I wanted to see what happens when we take a body with much more half-planes than dimensions, n << m.
Here with random direction walk.
Randomness has some cost. But the intersection (mostly dotProduct) is almost all of the work.

![1-vtune](./experiment1_cross/1_vtune.jpeg)
![1-pc](./experiment1_cross/1_pc.jpeg)
![1-vtune-2](./experiment1_cross/1_vtune_2.jpeg)


``` 
benchmark/benchmark_A1_volume -b "r=1,polytopeTranspose=false,generator=cross_rn_8" -c "step_size=1000000" -f "walk_f=walkCoord_ref"
```
Same with intersectCoord. Most of cost goes on intersectCoord. 
VTune measures the cost of cacheUpdateCoord surprisingly low, compared to pc calculations.
![2-pc](./experiment1_cross/2_pc.jpeg)
![2-vtune](./experiment1_cross/2_vtune.jpeg)
![2-vtune-2](./experiment1_cross/2_vtune_2.jpeg)

## Cube Rotated


``` 
benchmark/benchmark_A1_volume -b "r=1,polytopeTranspose=true,generator=cube_rot_r1.0_20" -c "step_size=1000000" -f "walk_f=walkCoord_ref"
``` 
Now on ballanced m = 2n. we see the logs at the end of each MCMC step is now visible (determine what shell the sample point falls into).
Again, vtune says cacheUpdateCoord is not so expensive...

![1-pc](./experiment2_cube/1_pc.jpeg)
![1-vtune](./experiment2_cube/1_vtune.jpeg)


``` 
-b "r=1,polytopeTranspose=true,generator=cube_rot_r1.0_100" -c "step_size=1000000" -f "walk_f=walkCoord_ref"
``` 
More dimensions, now with PolytopeT. intersectCoord is basically all of the work according to vTune.
Seems lots of time is spent on that one conditional... not sure if vtune measures this right.

![2-pc](./experiment2_cube/2_pc.jpeg)
![2-vtune](./experiment2_cube/2_vtune.jpeg)
![2-vtune-2](./experiment2_cube/2_vtune_2.jpeg)


``` 
benchmark/benchmark_A1_volume -b "r=1,polytopeTranspose=true,generator=cube_rot_r1.0_10" -c "step_size=1000000" -f "walk_f=walk_ref"
``` 
Back to random direction walk. 10-dim.
Cost of random is really expensive for random direction walk!

![3-pc](./experiment2_cube/3_pc.jpeg)
![3-vtune](./experiment2_cube/3_vtune.jpeg)
![3-vtune-2](./experiment2_cube/3_vtune_2.jpeg)

```
-b "r=1,polytopeTranspose=false,generator=cube_rot_r1.0_100" -c "step_size=10000" -f "walk_f=walk_ref" 
```
Also high costs for randomness with n=100, but somehow the intersect is now more expensive.
![4-pc](./experiment2_cube/4_pc.jpeg)
![4-vtune](./experiment2_cube/4_vtune.jpeg)




 
