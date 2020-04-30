# Experiments Log

## Cross

´´´
benchmark/benchmark\_A1\_volume -b "r=1,polytopeTranspose=false,generator=cross\_rn\_8" -c "step\_size=1000000"
´´´ 
![1-vtune](./experiment1_cross/1_vtune.jpeg)
![1-pc](./experiment1_cross/1_pc.jpeg)
![1-vtune-2](./experiment1_cross/1_vtune_2.jpeg)


´´´ 
benchmark/benchmark\_A1\_volume -b "r=1,polytopeTranspose=false,generator=cross\_rn\_8" -c "step\_size=1000000" -f "walk\_f=walkCoord\_ref"
´´´ 
VTune measures the cost of cacheUpdateCoord surprisingly low, compared to pc calculations.
![2-pc](./experiment1_cross/2_pc.jpeg)
![2-vtune](./experiment1_cross/2_vtune.jpeg)
![2-vtune-2](./experiment1_cross/2_vtune_2.jpeg)

## Cube Rotated


´´´ 
benchmark/benchmark\_A1\_volume -b "r=1,polytopeTranspose=true,generator=cube\_rot\_r1.0\_20" -c "step\_size=1000000" -f "walk\_f=walkCoord\_ref"
´´´ 

Again, vtune says cacheUpdateCoord is not so expensive...

![1-pc](./experiment2_cube/1_pc.jpeg)
![1-vtune](./experiment2_cube/1_vtune.jpeg)


´´´ 
-b "r=1,polytopeTranspose=true,generator=cube\_rot\_r1.0\_100" -c "step\_size=1000000" -f "walk\_f=walkCoord\_ref"
´´´ 
Seems lots of time is spent on that one conditional... not sure if vtune measures this right.

![2-pc](./experiment2_cube/2_pc.jpeg)
![2-vtune](./experiment2_cube/2_vtune.jpeg)
![2-vtune-2](./experiment2_cube/2_vtune_2.jpeg)


´´´ 
benchmark/benchmark\_A1\_volume -b "r=1,polytopeTranspose=true,generator=cube\_rot\_r1.0\_10" -c "step\_size=1000000" -f "walk\_f=walk\_ref"
´´´ 
Cost of random is really expensive for random direction walk!
![3-pc](./experiment2_cube/3_pc.jpeg)
![3-vtune](./experiment2_cube/3_vtune.jpeg)
![3-vtune-2](./experiment2_cube/3_vtune_2.jpeg)

´´´
-b "r=1,polytopeTranspose=false,generator=cube\_rot\_r1.0\_100" -c "step\_size=10000" -f "walk\_f=walk\_ref" 
´´´
Also high costs for randomness with n=100.
![4-pc](./experiment2_cube/4_pc.jpeg)
![4-vtune](./experiment2_cube/4_vtune.jpeg)




 
