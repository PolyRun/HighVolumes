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

typedef struct FTpair {
   double t0,t1;
} FTpair;

// This function must also reset/update the squaredNorm_cached function, which this intersectCoord depends on
typedef FTpair (*Ball_intersectCoord_cached_f_t)(const int n, const FT r, const FT* x, const int d, FT* cache);
extern Ball_intersectCoord_cached_f_t Ball_intersectCoord_cached;

FTpair Ball_intersectCoord_cached_ref(const int n, const FT r, const FT* x,const int d, FT* cache);

// ----------------------------- stuff for 4/8 way calculation:
typedef struct FTpair4 {
   __m256d low0,hi0;
} FTpair4;
typedef struct FTpair8 {
   __m256d low0,low1,hi0,hi1;
} FTpair8;
typedef struct FTset8 {
   __m256d set0,set1;
} FTset8;

__m256d squaredNorm_cached4(const FT* v, const int n, const FT* cache);
void squaredNorm_cached4_reset(const FT* v, const int n, FT* cache);
void squaredNorm_cached4_update(const FT* v, const int d, const __m256d dx, const int n, FT* cache);

FTset8 squaredNorm_cached8(const FT* v, const int n, const FT* cache);
void squaredNorm_cached8_reset(const FT* v, const int n, FT* cache);
void squaredNorm_cached8_update(const FT* v, const int d, const FTset8 dx, const int n, FT* cache);

FTpair4 Ball_intersectCoord_cached4(const int n, const FT r, const FT* x,const int d, FT* cache);
FTpair8 Ball_intersectCoord_cached8(const int n, const FT r, const FT* x,const int d, FT* cache);

// ------------------------------- arbitrary exponent numbers:

typedef struct ArbitraryExpNum {
   FT num; // regular number, may go to inf or zero if out of range
   FT numLog; // logarithm of number, should stay consistent, but loose a bit of precision.
} ArbitraryExpNum;

ArbitraryExpNum ArbitraryExpNum_new(FT val);
ArbitraryExpNum ArbitraryExpNum_mul(ArbitraryExpNum a, FT val);
ArbitraryExpNum ArbitraryExpNum_mul2(ArbitraryExpNum a, ArbitraryExpNum b);
void ArbitraryExpNum_print(ArbitraryExpNum a);




typedef void (*shell_cache_init_f_t)(FT *cache, FT r0, int l, FT stepFac);
extern shell_cache_init_f_t shell_cache_init;

void shell_cache_init_nocache(FT *cache, FT r0, int l, FT stepFac);
void shell_cache_init_ref(FT *cache, FT r0, int l, FT stepFac);

typedef int (*shell_idx_f_t)(FT x2, FT r0, FT stepFac, int l, FT *cache);
extern shell_idx_f_t shell_idx;

// default, doesn't use shell_cache -> does nothing
int shell_idx_nocache(FT x2, FT r0, FT stepFac, int l, FT *cache);
// simple linear search over the shells to avoid logs
int shell_idx_ref(FT x2, FT r0, FT stepFac, int l, FT *cache);
// binary search over shells
int shell_idx_binary(FT x2, FT r0, FT stepFac, int l, FT *cache);

#endif
