# Optimizations implemented

* dotProduct
  * 2 accumulators helps
  * input too small often - vectorization so far lead to no speedup
* Polytope vs PolytopeT
  * Column access speeds up intersectCoord (more locality)
  * expect it is easier to vectorize other functions
* intersectCoord(T): use cache to store dotProducts
  * details below

# Optimization Ideas

* dotProduct / vectorNorm
  * 2 accumulators brought some speedup
  * max 4 fps, fast around n=100, 2 fps below n=30
  * vectorization causes overhead to reduce final sum vector
  * could implement versions for different n?
  * vectorNorm should also be checked - not benchmarked yet

* Polytope.intersect
  * mostly dotProduct, could be fused to decrease reads.
  * may be harder to vectorize if dotProducts are subfunctions
  * conditionals / min / max could also be bottlenecks

* PolytopeT.intersect
  * vectorize. take k rows/constraints at a time, k multiple of 4
  * turn conditionals / min / max into vector instructions 

* Polytope(T).intersectCoord (with or without cache)
  * reduces flop/memory access -> about factor n
  * but: Polytope flop / memory density is now worse -> potential for improvement! 
  * and PolytopeT is much slower than Polytope for the impl without cache (dotProduct only fast for row-format)
  * but they are about the same for the cached version
  * So focus on improving PolytopeT.intersectCoord\_cached
![intersectCoord-cached](./optimizations/opt1_intersectCoord_cached_100.jpeg)
![intersectCoord-cached-T](./optimizations/opt1_intersectCoord_cached_100_T.jpeg)

* Parallelize to multiple walk points x
  * reduce intersection to MMM
  * probably produces lots of work (new signatures, more tests, etc)

