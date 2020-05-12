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

void PolytopeT_intersectCoord_vectorized(const void *p, const FT* x,
                                         const int d, FT *t0_out, FT *t1_out, void *cache);

#endif // HEADER_POLYTOPET_INTERSECTCOORD_H



