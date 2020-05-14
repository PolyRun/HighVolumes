#include "cacheUpdateCoord.h"

#include <immintrin.h>


void Ellipsoid_cacheUpdateCoord_c(const void* o, const int d, const FT dx, void* cache){
   Ellipsoid* e = (Ellipsoid*)o;
   const int n = e->n;
   FT* c = (FT*)cache;
   FT* a_p = Ellipsoid_get_a_p(e,d,d);

   __m128d cd_   = _mm_load_sd (&c[d]);
   __m128d two_ = _mm_set_sd (2.0);
   
   __m128d c2_  = _mm_mul_sd (two_, cd_);

   __m128d a_   = _mm_load_sd (a_p);
   __m128d dx_  = _mm_set_sd (dx);
   
   __m128d adx_ = _mm_mul_sd (a_, dx_);

   __m128d c2adx_ = _mm_add_sd (c2_, adx_);

   __m128d cn_   = _mm_load_sd (&c[n]);

   cn_ = _mm_fmadd_sd (dx_, c2adx_, cn_);

   _mm_store_sd (&c[n], cn_);

   for(int i=0; i<n; i++) {
      c[i] += dx * Ellipsoid_get_a(e,i,d);
   }
}

void Ellipsoid_cacheUpdateCoord_vec(const void* o, const int d, const FT dx, void* cache){
   Ellipsoid* e = (Ellipsoid*)o;
   const int n = e->n;
   FT* c = (FT*)cache;
   FT* a_p = Ellipsoid_get_a_p(e,d,d);

   __m128d cd_   = _mm_load_sd (&c[d]);
   __m128d two_ = _mm_set_sd (2.0);
   
   __m128d c2_  = _mm_mul_sd (two_, cd_);

   __m128d a_   = _mm_load_sd (a_p);
   __m128d dx_  = _mm_set_sd (dx);
   
   __m128d adx_ = _mm_mul_sd (a_, dx_);

   __m128d c2adx_ = _mm_add_sd (c2_, adx_);

   __m128d cn_   = _mm_load_sd (&c[n]);

   cn_ = _mm_fmadd_sd (dx_, c2adx_, cn_);

   _mm_store_sd (&c[n], cn_);

   FT* dcol = Ellipsoid_get_Ti(e, d);

   for(int i=0; i<n; i++) {
      c[i] += dx * Ellipsoid_get_a(e,i,d);
   }
}