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

void Ellipsoid_cacheUpdateCoord_fma(const void* o, const int d, const FT dx, void* cache){
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

   for(int i = 0; i < n; ++i) {
      __m128d dcoli_ = _mm_load_sd(&dcol[i]);
      __m128d ci_ = _mm_load_sd(&c[i]);
      ci_ = _mm_fmadd_sd(dx_, dcoli_, ci_);
      _mm_store_sd(&c[i], ci_);
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

   __m256d dxp_ = _mm256_set1_pd (dx);

   int i = 0;
   for(; i<n-3; i+=4) {
      __m256d dcoli_ = _mm256_load_pd(&dcol[i]);
      __m256d ci_ = _mm256_load_pd(&c[i]);
      ci_ = _mm256_fmadd_pd(dxp_, dcoli_, ci_);
      _mm256_store_pd(&c[i], ci_);
   }
   for(; i < n; ++i) {
      __m128d dcoli_ = _mm_load_sd(&dcol[i]);
      __m128d ci_ = _mm_load_sd(&c[i]);
      ci_ = _mm_fmadd_sd(dx_, dcoli_, ci_);
      _mm_store_sd(&c[i], ci_);
   }
}

void Ellipsoid_cacheUpdateCoord_vec_u2(const void* o, const int d, const FT dx, void* cache){
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

   FT* c2 = c+4;
   FT* dcol2 = dcol + 4;

   __m256d dxp_ = _mm256_set1_pd (dx);

   int i = 0;
   for(; i<n-7; i+=8) {
      __m256d dcoli_ = _mm256_load_pd(dcol);
      __m256d ci_ = _mm256_load_pd(c);
      ci_ = _mm256_fmadd_pd(dxp_, dcoli_, ci_);
      _mm256_store_pd(c, ci_);
      __m256d dcoli2_ = _mm256_load_pd(dcol2);
      __m256d ci2_ = _mm256_load_pd(c2);
      ci2_ = _mm256_fmadd_pd(dxp_, dcoli2_, ci2_);
      _mm256_store_pd(c2, ci2_);
      dcol += 8;
      dcol2 += 8;
      c+=8;
      c2+=8;
   }
   for(; i < n; ++i) {
      __m128d dcoli_ = _mm_load_sd(dcol);
      __m128d ci_ = _mm_load_sd(c);
      ci_ = _mm_fmadd_sd(dx_, dcoli_, ci_);
      _mm_store_sd(c, ci_);
      dcol++;
      c++;
   }
}

void Ellipsoid_cacheUpdateCoord_vec_u4(const void* o, const int d, const FT dx, void* cache){
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

   FT* c2 = c+4;
   FT* dcol2 = dcol + 4;
   FT* c3 = c+8;
   FT* dcol3 = dcol + 8;
   FT* c4 = c+12;
   FT* dcol4 = dcol + 12;

   __m256d dxp_ = _mm256_set1_pd (dx);

   int i = 0;
   for(; i<n-15; i+=16) {
      __m256d dcoli_ = _mm256_load_pd(dcol);
      __m256d ci_ = _mm256_load_pd(c);
      ci_ = _mm256_fmadd_pd(dxp_, dcoli_, ci_);
      _mm256_store_pd(c, ci_);
      __m256d dcoli2_ = _mm256_load_pd(dcol2);
      __m256d ci2_ = _mm256_load_pd(c2);
      ci2_ = _mm256_fmadd_pd(dxp_, dcoli2_, ci2_);
      _mm256_store_pd(c2, ci2_);
      __m256d dcoli3_ = _mm256_load_pd(dcol3);
      __m256d ci3_ = _mm256_load_pd(c3);
      ci3_ = _mm256_fmadd_pd(dxp_, dcoli3_, ci3_);
      _mm256_store_pd(c3, ci3_);
      __m256d dcoli4_ = _mm256_load_pd(dcol4);
      __m256d ci4_ = _mm256_load_pd(c4);
      ci4_ = _mm256_fmadd_pd(dxp_, dcoli4_, ci4_);
      _mm256_store_pd(c4, ci4_);
      dcol += 16;
      dcol2 += 16;
      c+=16;
      c2+=16;
   }
   for(; i < n; ++i) {
      __m128d dcoli_ = _mm_load_sd(dcol);
      __m128d ci_ = _mm_load_sd(c);
      ci_ = _mm_fmadd_sd(dx_, dcoli_, ci_);
      _mm_store_sd(c, ci_);
      dcol++;
      c++;
   }
}
