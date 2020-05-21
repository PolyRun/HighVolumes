#Analysis PolytopeJIT

Note: all of these measurements were taken on a machine that had some noise in the performance due to other applications running (eg browser). I tried to fix the situation, but we intend to run the benchmarks again on a different machine for the presentation and report.

## jit_test

We run ymm vectorized mul and max operation.
Max has gap of 1, thus we expect one every cycle, plus a mul.
4-way vectorized, so 8 flops per cycle at most.

Data test1 is linear data access, test2 random acces reads.
We do more reads for higher n, where also the data becomes larger.
We see growing performance and data io for higher n.

We can see that up to about n=100 the latency is dominating.
And we rise from about 5-20 bytes/cycle when going from n=3..100.

[Download jit_test io eps](https://gitlab.inf.ethz.ch/COURSE-ASL2020/team014/-/raw/master/optimizations/analysis_polytopeJIT/jit_benchmark/01_linear_vs_random_low_range/jit_test_io_mean.eps?inline=false)

[Download jit_test performance eps](https://gitlab.inf.ethz.ch/COURSE-ASL2020/team014/-/raw/master/optimizations/analysis_polytopeJIT/jit_benchmark/01_linear_vs_random_low_range/jit_test_performance_mean.eps?inline=false)

[Download jit_test runtime eps](https://gitlab.inf.ethz.ch/COURSE-ASL2020/team014/-/raw/master/optimizations/analysis_polytopeJIT/jit_benchmark/01_linear_vs_random_low_range/jit_test_runtime_mean.eps?inline=false)

## jit_test3 and jit_test4

Now we look at a case, where we generate a lot of little "kernels", little functions that do only a bit of work.
Here we expect high overhead for calling the functions.

jit_test3 is with ymm, jit_test4 with xmm registers. Thus we expect the latter to be half as fast, if the latency is dominant.

There are m kernels, and each kernel does w work (reads w ymm or xmm registers of data). The data is read from a random position in the n long array.

For ymm, we read about 0.5 bytes/cycle at least, and have 0.1 flops/cycle.
This increases both with w, but also if we have less kernels/functions.

This is probably due to instruction cache limits. Or jump prediction issues, when jumping to the right function.

[Download jit_test3 io eps](https://gitlab.inf.ethz.ch/COURSE-ASL2020/team014/-/raw/master/optimizations/analysis_polytopeJIT/jit_benchmark/02_ymm_kernels/jit_test_3_io_mean.eps?inline=false)

[Download jit_test3 performance eps](https://gitlab.inf.ethz.ch/COURSE-ASL2020/team014/-/raw/master/optimizations/analysis_polytopeJIT/jit_benchmark/02_ymm_kernels/jit_test_3_performance_mean.eps?inline=false)

[Download jit_test3 runtime eps](https://gitlab.inf.ethz.ch/COURSE-ASL2020/team014/-/raw/master/optimizations/analysis_polytopeJIT/jit_benchmark/02_ymm_kernels/jit_test_3_runtime_mean.eps?inline=false)

[Download jit_test4 io eps](https://gitlab.inf.ethz.ch/COURSE-ASL2020/team014/-/raw/master/optimizations/analysis_polytopeJIT/jit_benchmark/02_ymm_kernels/jit_test_4_io_mean.eps?inline=false)

[Download jit_test4 performance eps](https://gitlab.inf.ethz.ch/COURSE-ASL2020/team014/-/raw/master/optimizations/analysis_polytopeJIT/jit_benchmark/02_ymm_kernels/jit_test_4_performance_mean.eps?inline=false)

[Download jit_test4 runtime eps](https://gitlab.inf.ethz.ch/COURSE-ASL2020/team014/-/raw/master/optimizations/analysis_polytopeJIT/jit_benchmark/02_ymm_kernels/jit_test_4_runtime_mean.eps?inline=false)

Note: this case is supposed to model our sparse body access pattern.

## Performance for Dense bodies

We take an n-dimensional hypercube, and rotate it in space to get a dense body (the constraint matrix is dense).

### intersectCoord

For each position, we read Ainv and the cache entry.
We mul and max/min. Thus, we have 2 flops and 16 bytes data.

For single value access (single=1 double), we have the rax or data way. Rax means it was read from the instruction, as an 64 bit integer payload and loaded through rax to xmm registers (there is no direct way to load instruction payloads into xmm/ymm registers). The Rax option turns out to be slower. It bloats the instruction bytes, and does not utilize the separation of instruction and data cache.

The acc option is to use multiple accumulators. If we only have one for min and one for max, we can only use 2 out of three slots for max/min (3 lat, 1 gap). This option did not really lead to speedup. Maybe the data layout does not allow it, or maybe the instructions come in too slow to be handled more efficiently, we are not sure.

As soon as we use xmm or ymm registers, we use blocking. This means we allow one register only to hold negative or positive values, potentially having to run a block twice, once for positive values with min, once with negatives and max.

Using xmm or even ymm registers gives some speedup, though not exactly linear. This is probably due to the fact that we cannot guarantee that adjacent values have the same sign, but it happens sometimes, so some speedup.

We also tried reordering the constraints to have better locality in the data pattern, so rearrange to have more values with the same sign close to each other in a column. This is a best effort algorithm where we do not know how good it really is. Unfortunately, we cannot see a signifficant difference.

[Download io eps](https://gitlab.inf.ethz.ch/COURSE-ASL2020/team014/-/raw/master/optimizations/analysis_polytopeJIT/dense_cubeRot/dense_polytopeJIT_intersect_io_mean.eps?inline=false)

[Download performance eps](https://gitlab.inf.ethz.ch/COURSE-ASL2020/team014/-/raw/master/optimizations/analysis_polytopeJIT/dense_cubeRot/dense_polytopeJIT_intersect_performance_mean.eps?inline=false)

[Download runtime eps](https://gitlab.inf.ethz.ch/COURSE-ASL2020/team014/-/raw/master/optimizations/analysis_polytopeJIT/dense_cubeRot/dense_polytopeJIT_intersect_runtime_mean.eps?inline=false)

Explanation for Performance Roof:

* As n is small, the latency dominates.
* For large n, throughput dominates.
  * This is then limited by the instruction cache and op decoding/scheduling (not sure how much).
  * From the measurements in jit_test, we know that we most likely will not get over 20 bytes/cycle for n=100 or smaller.   
  * Further, some blocks have to be done twice, as they contain both positive and negative values. Hence, we might expect about 10 bytes/cycle.
  * We also get to about 1 flop/cycle, where according to jit_test we could have at most a bit over 2 flops/cycle.
  * This has to be compared to the theoretical 8 flops/cycle the instruction mix of mul and max/min would allow.

### cacheUpdateCoord

Here, we read/write the cache, updating the dot-product in one dimension. This is a simple fmadd, and 3 doubles io.

This operation is still memory bound, though now we have removed the restriction of min/max gaps as in the intersectCoord. Thus, we also only need to look at each block exactly once.

This should make it signifficantly simpler from an op scheduling perspective.

Thus we have two factors to the faster performance for this function: simpler op scheduling and handling each block exactly once.

[Download io eps](https://gitlab.inf.ethz.ch/COURSE-ASL2020/team014/-/raw/master/optimizations/analysis_polytopeJIT/dense_cubeRot/dense_polytopeJIT_update_io_mean.eps?inline=false)

[Download performance eps](https://gitlab.inf.ethz.ch/COURSE-ASL2020/team014/-/raw/master/optimizations/analysis_polytopeJIT/dense_cubeRot/dense_polytopeJIT_update_performance_mean.eps?inline=false)

[Download runtime eps](https://gitlab.inf.ethz.ch/COURSE-ASL2020/team014/-/raw/master/optimizations/analysis_polytopeJIT/dense_cubeRot/dense_polytopeJIT_update_runtime_mean.eps?inline=false)

## Performance for Sparse Bodies


