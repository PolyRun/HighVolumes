#ifndef HEADER_POLYTOPECSC_INTERSECTCOORD_H
#define HEADER_POLYTOPECSC_INTERSECTCOORD_H

#include "../polytopeCSC.h"

void PolytopeCSC_intersectCoord_cached_ref(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache);
void PolytopeCSC_intersectCoord_cached_withb(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache);
void PolytopeCSC_intersectCoord_cached_vec(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache);
void PolytopeCSC_intersectCoord_cached_vec_inline(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache);


#endif // HEADER_POLYTOPECSC_INTERSECTCOORD_H



