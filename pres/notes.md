# Presentation Notes
* quickly explain title

## 1
* Convex body, n-dimensional, want to get volume

## 2
* estimate PI, sample uniformly in square
* check what falls into circle
* get ratio of volumes, calculate area of circle
* New Problem: curse of dimensionality
* Need exponential #points to hit ball

## 3
* For this talk: consider only Polytopes
  * Set of linear constraints: Ax <= b
  * Here: 5 linear inequalities
* Subdivision: solution to "New Problem", curse of dimensionality
* Assume: have inner and outer ball given (explain what contains what)
  * there exists algo, out of scope
* Zone is intersection

## 4
* solved curse of dimensionality
* ratio at most 2 -> probability to hit smaller zone at least 1/2, not exp.
* find ratio: via random sample points

## 5
* Telescoping product for volume
* ratios, volume Zone1, cancel, Zone_n == Body

## 6
* random COORD direction, intersection, random point on segment
* repeat many times
* approximation for uniform random point in Zone
* Bulk of work is computing intersection, rest is comparatively cheap.

## 7 [exercize this well !]
* How to find segment: compute all intersections
* Get distance to each constraint, find closest in both dir

## 8
* Look at computatin of t_i, need to optimize
* 2 dot products -> EXPENSIVE !
* lower is single entry, d is unit vec in coord dir
* Update x: [red] one entry changes
* no need to recompute dot prod: cache and subtract
* no dot product left

## 9
* Implemented, optimized and benchmarked.
* Haswell, customary compiler flags, signifficant L1 cache

## 10
* Most work is intersection computation
* Here: see benefit of cacheing
* Cube, various dim, rotated, intersection 2*dim constraints
* Measured Runtime
* Error bars: 95% Confidence Interval (top,bottom 2.5% discarded)

## 11
* Zoom in: top-to-bottom, speedup
* div slow, precompute for each entry in constraint Matrix, just mul

## 12
* Performance: flops/cycle, bottum-up
* cached_vec: div performance bound below 1 flops/cycle
* precomp inv: reduced to mul, bound at 5 flops/cycle
* Why only 3 reached? Cube does not have big enough columns.
* Cost of latency. More constraints gets better.

## 13
* Interesting special case.
* Example: only 2 variables per constraint.
* Sparse constraint matrix.
* How to implement?

## 14
* Sparse Column format, as in lecture.
* Non-zeros contiguous, plus row indices.
* Note: relevant entries in cache NOT contiguous.
* Store inv, just load and mul.

* Turns out: path depends only on A.
* Repeat intersection many times.
* can precompute branching, left/right?

## 15
* Yes!
* Second approach: kernel/little funcion.
* for each coord dir: have an intersection kernel, compute intersection.
* code at runtime, depending on A.

## 16
* Want to see how good for different sparsity levels.
* Runtime intersect for 100-dim, 1000 constraints.
* X-axis: % non-zero elem in A, log-scale.
* Blue line: our best dense impl, as discussed earlier.
* 2 classes of sparse impl, CSC and JIT.
* top-down: optimizations for CSC, data access for JIT.
* dense stays constant, sparse only good for sparse enough data.


## 17
* Roofline for the same experiment.
* Roofs: data and performance
* Name: JIT, CSC, dense_impl
  * sparse-to-dense: dense -> more ops -> lower impact of latency
* CSC red: div bound
* dense_impl: much more work

## 18
* Runtime, for 2-var, 10*dim constraints.
* optimized CSC about as good as JIT, but in very sparse JIT can be better.
* CSC must also load indices, and more ops.

## 19
* What we learned:
* Cache to reduce work.
* Use sparsity. Many reps (temp locality) -> kernel idea. Limit: L1i cache.




