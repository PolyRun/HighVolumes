
#ifndef HEADER_VOLUME_COST_HPP
#define HEADER_VOLUME_COST_HPP

#include "../util/performance_counter.hpp"

typedef void (*xyz_cost_f) (const int);
void xyz_f1_cost(const int n);

// ----------------------------------------------- random
typedef void (*random_int_cost_f)(const void* o);
typedef void (*random_int_in_range_cost_f)(const void* o);
typedef void (*random_double_in_range_cost_f)(const void* o);
typedef void (*random_double_0_1_cost_f)(const void* o);
typedef void (*random_double_normal_cost_f)(const void* o);
void Random_int_cost_ref(const void * o);
void Random_int_in_range_cost_ref(const void * o);
void Random_double_in_range_cost_ref(const void * o);
void Random_double_0_1_cost_ref(const void * o);
void Random_double_normal_cost_ref(const void * o);

typedef void (*rand256d_cost_f_t)();
void sr_rand256d_cost_ref();

// ----------------------------------------------- linalg
typedef void (*dotProduct_cost_f)(const int);
void dotProduct_cost_ref(const int n);
typedef void (*squaredNorm_cost_f)(const int);
void squaredNorm_cost_ref(const int n);

typedef void (*Ball_intersectCoord_cost_f)(const int);
void Ball_intersectCoord_cost_ref(const int n);
typedef void (*Ball_intersect_cost_f)(const int);
void Ball_intersect_cost_ref(const int n);

typedef void (*Ball_intersectCoord_cached_cost_f)(const int);
void Ball_intersectCoord_cached_cost_ref(const int n);
void Ball_intersectCoord_cached4_cost_ref(const int n);
void Ball_intersectCoord_cached8_cost_ref(const int n);

// ----------------------------------------------- Body
typedef void (*intersect_cost_f)(const void*);
typedef void (*intersectCoord_cost_f)(const void*);
typedef void (*cacheUpdateCoord_cost_f)(const void*);
typedef void (*cacheReset_cost_f)(const void*);

// ----------------------------------------------- Polytope
void Polytope_intersect_cost_ref(const void* o);
void Polytope_intersectCoord_cost_ref(const void* o);
void Polytope_intersectCoord_cached_cost_ref(const void* o);
void Polytope_cacheUpdateCoord_cost_ref(const void* o);
void Polytope_cacheReset_cost_ref(const void* o);

// ----------------------------------------------- PolytopeT
void PolytopeT_intersect_cost_ref(const void* o);
void PolytopeT_intersectCoord_cost_ref(const void* o);
void PolytopeT_intersectCoord_cost_vec(const void* o);
void PolytopeT_intersectCoord_cached_cost_ref(const void* o);
void PolytopeT_intersectCoord_cached_b_cost_ref(const void* o);
void PolytopeT_cacheUpdateCoord_cost_ref(const void* o);
void PolytopeT_cacheUpdateCoord_b_cost_ref(const void* o);
void PolytopeT_cacheUpdateCoord_b_cost_vec(const void* o);
void PolytopeT_cacheReset_cost_ref(const void* o);
void PolytopeT_cacheReset_b_cost_ref(const void* o);
void PolytopeT_cacheReset_b_cost_vec(const void* o);

void PolytopeT_intersectCoord4_cost_ref(const void* o);
void PolytopeT_intersectCoord8_cost_ref(const void* o);
void PolytopeT_cacheUpdateCoord4_cost_ref(const void* o);
void PolytopeT_cacheUpdateCoord8_cost_ref(const void* o);
void PolytopeT_cacheReset4_cost_ref(const void* o);
void PolytopeT_cacheReset8_cost_ref(const void* o);


// ----------------------------------------------- Ellipsoid
void Ellipsoid_intersect_cost_ref(const void* o);
void Ellipsoid_intersectCoord_cost_ref(const void* o);
void Ellipsoid_intersectCoord_cached_cost_ref(const void* o);
void Ellipsoid_intersectCoord_cached_cost_reord_fma(const void* o);
void Ellipsoid_cacheUpdateCoord_cost_ref(const void* o);
void Ellipsoid_cacheReset_cost_ref(const void* o);

// ----------------------------------------------- PolytopeCSC


typedef void (*mvm_cost_f)(const PolytopeCSC *p);

// NOTE: the actual #flops & #bytes depends on direction d (c.f. nonzerosCSC)
// we average #non-zeros over all columns to get estimate on cost of these functions
void PolytopeCSC_intersect_cost_ref(const void *o);
void PolytopeCSC_mvm_cost(const PolytopeCSC *p);
void PolytopeCSC_intersectCoord_cost_ref(const void *o);
void PolytopeCSC_intersectCoord_cached_cost_ref(const void *o);
void PolytopeCSC_cacheReset_cost_ref(const void *o);
void PolytopeCSC_cacheUpdateCoord_cost_ref(const void *o);
void PolytopeCSC_intersectCoord_cached_cost_withb(const void *o);
void PolytopeCSC_intersectCoord_cached_cost_vec(const void *o);
void PolytopeCSC_cacheReset_cost_withb(const void *o);
void PolytopeCSC_cacheUpdateCoord_cost_withb(const void *o);

void PolytopeCSC_intersectCoord4_cost_ref(const void* o);
void PolytopeCSC_intersectCoord8_cost_ref(const void* o);
void PolytopeCSC_cacheUpdateCoord4_cost_ref(const void* o);
void PolytopeCSC_cacheUpdateCoord8_cost_ref(const void* o);
void PolytopeCSC_cacheReset4_cost_ref(const void* o);
void PolytopeCSC_cacheReset8_cost_ref(const void* o);


// ----------------------------------------------- PolytopeJIT
void PolytopeJIT_intersect_cost_ref(const void* o);
void PolytopeJIT_intersectCoord_cost_ref(const void* o);
void PolytopeJIT_cacheUpdateCoord_cost_ref(const void* o);
void PolytopeJIT_cacheReset_cost_ref(const void* o);

void PolytopeJIT_intersectCoord4_cost_ref(const void* o);
void PolytopeJIT_intersectCoord8_cost_ref(const void* o);
void PolytopeJIT_cacheUpdateCoord4_cost_ref(const void* o);
void PolytopeJIT_cacheUpdateCoord8_cost_ref(const void* o);
void PolytopeJIT_cacheReset4_cost_ref(const void* o);
void PolytopeJIT_cacheReset8_cost_ref(const void* o);

// ----------------------------------------------- volume

typedef void (*volume_cost_f)(const int, const int, const void**,const Body_T**);
void volume_cost_ref(const int n, const int bcount, const void** body, const Body_T** type);
void volume_coord_1_cost_ref(const int n, const int bcount, const void** body, const Body_T** type);
void volume_coord_4_cost_ref(const int n, const int bcount, const void** body, const Body_T** type);
void volume_coord_8_cost_ref(const int n, const int bcount, const void** body, const Body_T** type);

typedef void (*walk_cost_f)(const int, int bcount, const void**, const Body_T**);
void walk_cost_ref(const int n, int bcount, const void** body, const Body_T** type);
void walkCoord_cost_ref(const int n, int bcount, const void** body, const Body_T** type);
void walkCoord_coord_1_cost_ref(const int n, int bcount, const void** body, const Body_T** type);
void walkCoord_coord_4_cost_ref(const int n, int bcount, const void** body, const Body_T** type);
void walkCoord_coord_8_cost_ref(const int n, int bcount, const void** body, const Body_T** type);


#endif // HEADER_VOLUME_COST_HPP
