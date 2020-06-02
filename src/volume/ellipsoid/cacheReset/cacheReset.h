#ifndef HEADER_ELLIPSOID_CACHERESET_H
#define HEADER_ELLIPSOID_CACHERESET_H

#include "../ellipsoid.h"

void Ellipsoid_cacheReset_reord(const void* o, const FT* x, void* cache);

void Ellipsoid_cacheReset_fma(const void* o, const FT* x, void* cache);

//void Ellipsoid_cacheReset_vec(const void* o, const FT* x, void* cache);

#endif // HEADER_ELLIPSOID_CACHERESET_H
