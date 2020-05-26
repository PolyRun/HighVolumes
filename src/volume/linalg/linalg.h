#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <immintrin.h>
#include <assert.h>

#include "../../random/prng.h"
#include "cholesky.h"

#ifndef LINALG_H
#define LINALG_H

// dotProduct:
//   assume no memory allignment!
typedef FT (*dotProduct_f_t)(const FT* u, const FT* v, const int n);
extern dotProduct_f_t dotProduct;
typedef FT (*squaredNorm_f_t)(const FT* v, const int n);
extern squaredNorm_f_t squaredNorm;
#include "dotProduct/dotProduct.h"

// intersect line x + t*d
//           with ball(0,r)
// assume x is in ball
// return t0,t1 for intersections
typedef void (*Ball_intersect_f_t)(const int n, const FT r, const FT* x, const FT* d, FT* t0, FT* t1);
extern Ball_intersect_f_t Ball_intersect;
void Ball_intersect_ref(const int n, const FT r, const FT* x, const FT* d, FT* t0, FT* t1);

typedef void (*Ball_intersectCoord_f_t)(const int, const FT, const FT*, const int, FT*, FT*);
extern Ball_intersectCoord_f_t Ball_intersectCoord;
void Ball_intersectCoord_ref(const int n, const FT r, const FT* x, const int d, FT* t0, FT* t1);

    
// calculate volume exactly for n-dim ball with radius r
FT Ball_volume(const int n, const FT r);

// given n elements of b bytes, want to get smallest multiple of 32 bytes that fits this.
int ceil_cache(const int n, const int b);

// --------------- cached squaredNorm:
typedef FT (*squaredNorm_cached_f_t)(const FT* v, const int n, const FT* cache);
extern squaredNorm_cached_f_t squaredNorm_cached;
typedef void (*squaredNorm_cached_reset_f_t)(const FT* v, const int n, FT* cache);
extern squaredNorm_cached_reset_f_t squaredNorm_cached_reset;
typedef void (*squaredNorm_cached_update_f_t)(const FT* v, const int d, const FT dx, const int n, FT* cache);
extern squaredNorm_cached_update_f_t squaredNorm_cached_update;
// Update is called after v[d] += dx;

FT squaredNorm_cached_ref(const FT* v, const int n, const FT* cache);
void squaredNorm_cached_reset_ref(const FT* v, const int n, FT* cache);
void squaredNorm_cached_update_ref(const FT* v, const int d, const FT dx, const int n, FT* cache);

FT squaredNorm_cached_refc(const FT* v, const int n, const FT* cache);
void squaredNorm_cached_reset_refc(const FT* v, const int n, FT* cache);
void squaredNorm_cached_update_refc(const FT* v, const int d, const FT dx, const int n, FT* cache);



#endif
