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

TODO
