#ifndef HEADER_POLYTOPET_INTERSECTCOORD_H
#define HEADER_POLYTOPET_INTERSECTCOORD_H

#include "../polytopeT.h"

void PolytopeT_intersectCoord_cached_ref(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache);

// nc: no conditional
// the if/continue statement seems to be very costly
// let's remove it
//
// nc1: forked from cached_ref
void PolytopeT_intersectCoord_cached_nc1(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache);
void PolytopeT_intersectCoord_cached_b_ref(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache);
void PolytopeT_intersectCoord_cached_b_vec(const void *p, const FT *x, const int d,FT *t0_out, FT *t1_out, void *cache);
void PolytopeT_intersectCoord_cached_b_vec2(const void *p, const FT *x, const int d,FT *t0_out, FT *t1_out, void *cache);
void PolytopeT_intersectCoord_cached_b_vec_inl(const void *p, const FT *x, const int d,FT *t0_out, FT *t1_out, void *cache);
void PolytopeT_intersectCoord_vectorized(const void *p, const FT* x, const int d, FT *t0_out, FT *t1_out, void *cache);

void PolytopeT_intersectCoord_cached_b_inv_ref(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache);
void PolytopeT_intersectCoord_cached_b_inv_vec(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache);
void PolytopeT_intersectCoord_cached_b_inv_vec_inl(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache);

void PolytopeT_cacheReset_b_ref(const void* o, const FT* x, void* cache);
void PolytopeT_cacheReset_b_vec(const void *p, const FT *x, void *cache);

void PolytopeT_cacheUpdateCoord_b_ref(const void* o, const int d, const FT dx, void* cache);
void PolytopeT_cacheUpdateCoord_b_vec(const void *p, const int d, const FT dx, void *cache);



FTpair4 PolytopeT_intersectCoord4_ref(const void* o, const FT* x, const int d, void* cache);
FTpair8 PolytopeT_intersectCoord8_ref(const void* o, const FT* x, const int d, void* cache);
void PolytopeT_cacheReset4_ref(const void* o, const FT* x, void* cache);
void PolytopeT_cacheReset8_ref(const void *p, const FT *x, void *cache);
void PolytopeT_cacheUpdateCoord4_ref(const void* o, const int d, const __m256d dx, void* cache);
void PolytopeT_cacheUpdateCoord8_ref(const void* o, const int d, const FTset8 dx, void* cache);

#endif // HEADER_POLYTOPET_INTERSECTCOORD_H



