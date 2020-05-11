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

   const FT aInv = 0.5/a;     // Issued at 0, available at 14

   FT b = 2.0*Az_c[d];        // Issued at 0, available at 4
   FT c = 4.0*Az_c[n];        // Issued at 0, available at 4

   FT bsquared = b*b;         // Issued at 4, available at 8
   
   // find t:
   const FT discr = bsquared - a*c; // Issued at 8(12), available at 16(12),    fNmadd here -4
   assert(discr >= 0);
   const FT sqrtDiscr = sqrt(discr); //Issued at 16, available at 34

   FT tmp = -b * aInv; // Issued at 14, available at 18

   *t0 = tmp - sqrtDiscr * aInv; // Issued at 34(38), available at 42(38),     fNmadd here -4
   *t1 = tmp + sqrtDiscr * aInv; // Issued at 34(38), available at 42(38),     fmadd here  -4
}

void Ellipsoid_intersectCoord_cached_reord_fma(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache) {
   Ellipsoid* e = (Ellipsoid*)o;
   const int n = e->n;
   FT* Az_c = (FT*)cache;

   FT* a_p = Ellipsoid_get_a_p(e,d,d);

   __m128d a = _mm_load_sd (a_p);
   __m128d half = _mm_set_sd (0.5);

   const __m128d aInv = _mm_div_sd (a, half);

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
