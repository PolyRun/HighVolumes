#Analysis PolytopeJIT

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

## Performance for Dense bodies:


