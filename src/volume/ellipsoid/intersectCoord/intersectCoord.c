#include "intersectCoord.h"

#include <immintrin.h>
#include <emmintrin.h>


void Ellipsoid_intersectCoord_cached_reord(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache) {
   Ellipsoid* e = (Ellipsoid*)o;
   const int n = e->n;
   FT* Az_c = (FT*)cache;

   FT a = Ellipsoid_get_a(e,d,d);

   const FT aInv = 0.5/a;

   FT b = 2.0*Az_c[d];
   FT c = Az_c[n];
   
   // find t:
   const FT discr = b*b - 4.0*a*c;
   assert(discr >= 0);
   const FT sqrtDiscr = sqrt(discr);

   *t0 = (-b - sqrtDiscr) * aInv;
   *t1 = (-b + sqrtDiscr) * aInv;
}

void Ellipsoid_intersectCoord_cached_reord2(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache) {
   Ellipsoid* e = (Ellipsoid*)o;
   const int n = e->n;
   FT* Az_c = (FT*)cache;

   FT a = Ellipsoid_get_a(e,d,d);

   const FT aInv = 0.5/a;

   FT b = 2.0*Az_c[d];
   FT c = 4.0*Az_c[n];

   FT bsquared = b*b;
   
   // find t:
   const FT discr = bsquared - a*c;
   assert(discr >= 0);
   const FT sqrtDiscr = sqrt(discr);

   *t0 = (-b - sqrtDiscr) * aInv;
   *t1 = (-b + sqrtDiscr) * aInv;
   
}

void Ellipsoid_intersectCoord_cached_reord3(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache) {
   Ellipsoid* e = (Ellipsoid*)o;
   const int n = e->n;
   FT* Az_c = (FT*)cache;

   FT a = Ellipsoid_get_a(e,d,d);

   const FT aInv = 0.5/a;

   FT b = 2.0*Az_c[d];
   FT c = 4.0*Az_c[n];

   FT bsquared = b*b;
   
   // find t:
   const FT discr = bsquared - a*c;
   assert(discr >= 0);
   const FT sqrtDiscr = sqrt(discr);

   FT tmp = -b * aInv;

   *t0 = tmp - sqrtDiscr * aInv;
   *t1 = tmp + sqrtDiscr * aInv;
}

/**
 * Reasoning:
 *                                   Skylake               Haswell
 * resVar   | Op             | IssuedAt | Finished | IssuedAt | Finished
 * ---------------------------------------------------------------------
 * aInv     | 0.5/a          |     0    |    14    |     0    |  14-20
 * b        | 2.0*Az_c[d]    |     0    |     4    |     0    |      5
 * c        | 4.0*Az_c[n]    |     0    |     4    |     0    |      5
 * bsquared | b*b            |     4    |     8    |     5    |     10
 * discr    | bsq - a*c      |     8    |    12    |    10    |     15
 * sqrtD    | sqrt(discr)    |    12    |    30    |    15    |     35
 * tmp      | -b*aInv        |    14    |    18    |  14-20   |  19-25
 * *t0      | tmp-sqrtD*aInv |    30    |    34    |    35    |     40
 * *t1      | tmp+sqrtD*aInv |    30    |    34    |    35    |     40
 * 
 * 11 flops in 34(40) cycles
 * do not count sqrtD*aInv twice!
 *
 * with cacheUpdate:
 * n=100:
 * 217 flops in 150 cycles -> 1.44 flops per cycle (measured: 1.65 - would be 248 flops or 131 cycles)
 **/

void Ellipsoid_intersectCoord_cached_reord_fma(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache) {
   Ellipsoid* e = (Ellipsoid*)o;
   const int n = e->n;
   FT* Az_c = (FT*)cache;

   FT* a_p = Ellipsoid_get_a_p(e,d,d);

   __m128d a = _mm_load_sd (a_p);
   __m128d half = _mm_set_sd (0.5);

   const __m128d aInv = _mm_div_sd (half, a);

   __m128d b = _mm_load_sd (&Az_c[d]);
   __m128d two = _mm_set_sd (2.0);
   b = _mm_mul_sd (two, b);

   __m128d c = _mm_load_sd (&Az_c[n]);
   __m128d four = _mm_set_sd (4.0);
   c = _mm_mul_sd (four, c);

   __m128d bsquared = _mm_mul_sd (b, b);

   // find t:
   const __m128d discr = _mm_fnmadd_sd (a, c, bsquared);
   __m128d zero = _mm_set_sd (0.0);
   const __m128d sqrtDiscr = _mm_sqrt_sd (zero, discr);

   __m128d tmp = _mm_fnmadd_sd (b, aInv, zero);

   __m128d r0 = _mm_fnmadd_sd (sqrtDiscr, aInv, tmp);
   __m128d r1 = _mm_fmadd_sd (sqrtDiscr, aInv, tmp);

   _mm_store_sd (t0, r0);
   _mm_store_sd (t1, r1);
}
