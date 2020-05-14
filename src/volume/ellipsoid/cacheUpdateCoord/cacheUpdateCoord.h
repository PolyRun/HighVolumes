#ifndef HEADER_ELLIPSOID_CACHEUPDATECOORD_H
#define HEADER_ELLIPSOID_CACHEUPDATECOORD_H

#include "../ellipsoid.h"

void Ellipsoid_cacheUpdateCoord_c(const void* o, const int d, const FT dx, void* cache);

void Ellipsoid_cacheUpdateCoord_fma(const void* o, const int d, const FT dx, void* cache);

void Ellipsoid_cacheUpdateCoord_vec(const void* o, const int d, const FT dx, void* cache);

void Ellipsoid_cacheUpdateCoord_vec_u2(const void* o, const int d, const FT dx, void* cache);

void Ellipsoid_cacheUpdateCoord_vec_u4(const void* o, const int d, const FT dx, void* cache);

#endif // HEADER_ELLIPSOID_CACHEUPDATECOORD_H
