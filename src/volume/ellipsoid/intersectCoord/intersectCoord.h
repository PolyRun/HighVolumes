#ifndef HEADER_ELLIPSOID_INTERSECTCOORD_H
#define HEADER_ELLIPSOID_INTERSECTCOORD_H

#include "../ellipsoid.h"

// Reordered instuctions
void Ellipsoid_intersectCoord_cached_reord(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache);

void Ellipsoid_intersectCoord_cached_reord2(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache);

void Ellipsoid_intersectCoord_cached_reord3(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache);

// Reordered instructions using FMAs
void Ellipsoid_intersectCoord_cached_reord_fma(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache);

#endif // HEADER_ELLIPSOID_INTERSECTCOORD_H
