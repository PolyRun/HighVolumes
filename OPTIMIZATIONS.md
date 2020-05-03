# Optimizations

* dotProduct / vectorNorm
  * 2 accumulators brought some speedup
  * input too small often - vectorization so far lead to no speedup
  * max 4 fps, fast around n=100, 2 fps below n=30
  * vectorization causes overhead to reduce final sum vector
  * could implement versions for different n?
  * vectorNorm should also be checked - not benchmarked yet

* Polytope vs PolytopeT
  * row vs column matrix A. allows for different impl, especially when vectorizing

* Polytope.intersect
  * mostly dotProduct, could be fused to decrease reads.
  * may be harder to vectorize if dotProducts are subfunctions
  * conditionals / min / max could also be bottlenecks

* PolytopeT.intersect
  * vectorize. take k rows/constraints at a time, k multiple of 4
  * turn conditionals / min / max into vector instructions 

* Polytope(T).intersectCoord (with or without cache - store dot product for polytope intersection)
  * reduces flop/memory access -> about factor n
  * but: Polytope flop / memory density is now worse -> potential for improvement! 
  * and PolytopeT is much slower than Polytope for the impl without cache (dotProduct only fast for row-format)
  * but they are about the same for the cached version
  * So focus on improving PolytopeT.intersectCoord\_cached
  * do not forget the costs of cacheUpdateCoord
    * for Polytope, this may be bad bc strided memory access, in PolytopeT this is a simple vector += vector x scalar (vfma)
![intersectCoord-cached](./optimizations/opt1_intersectCoord_cached_100.jpeg)
![intersectCoord-cached-T](./optimizations/opt1_intersectCoord_cached_100_T.jpeg)

* Ellipsoid
  * Cost difference is not so big between intersect and intersectCoord, only flop count is halved.
  * probably could gain some speedup with traditional MVM techniques
  * Maybe there could be a way to cache the MVM? Maybe in the intersectCoord case this could lead to something.
![ellipsoid-intersect](./optimizations/opt1_intersect_ellipsoid_100.jpeg)


* Parallelize to multiple walk points x
  * reduce intersection to MMM
  * probably produces lots of work (new signatures, more tests, etc)

