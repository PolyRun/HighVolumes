#include "cacheReset.h"

#include <immintrin.h>


void Ellipsoid_cacheReset_reord(const void* o, const FT* x, void* cache){
   const Ellipsoid* e = (Ellipsoid*)o;
   FT* c = (FT*)cache;
   const int n = e->n;

   // compute x[index]-e->a[index]
   FT xea[n];
   for(int index = 0; index < n; ++index) {
      xea[index] = x[index] - e->a[index];
   }

   c[n] = -1.0;
   for(int i=0;i<n;i++) {
      const FT* Ai = Ellipsoid_get_Ai(e,i);
      c[i] = 0;
      for(int j=0; j<n; j++) {
         c[i] += Ai[j] * xea[j];
      }
      c[n] += xea[i] * c[i];
   }
}

void Ellipsoid_cacheReset_fma(const void* o, const FT* x, void* cache){
   const Ellipsoid* e = (Ellipsoid*)o;
   FT* c = (FT*)cache;
   const int n = e->n;

   // compute x[index]-e->a[index]
   __m128d xea[n];
   for(int index = 0; index < n; ++index) {
      __m128d x_index  = _mm_load_sd(&x[index]);
      __m128d ea_index = _mm_load_sd(&e->a[index]);
      xea[index] = _mm_add_sd(x_index, ea_index);
   }

   __m128d cn = _mm_set_sd(-1.0);
   for(int i=0;i<n;i++) {
      const FT* Ai = Ellipsoid_get_Ai(e,i);
      __m128d ci = _mm_set_sd(0.0);
      for(int j=0; j<n; j++) {
        __m128d Aij = _mm_load_sd(&Ai[j]);
        ci = _mm_fmadd_sd(Aij, xea[j], ci);
      }
      _mm_store_sd(&c[i], ci);
      cn = _mm_fmadd_sd(ci, xea[i], cn);
   } 
   _mm_store_sd(&c[n], cn);
}


/*
   Not considered further since it uses interleaved SSE and AVX together with reduction on __256ds
*/
/*void Ellipsoid_cacheReset_vec(const void* o, const FT* x, void* cache){
   const Ellipsoid* e = (Ellipsoid*)o;
   FT* c = (FT*)cache;
   const int n = e->n;

   // compute x[index]-e->a[index]
   __m256d xea[n>>2];
   int index = 0;
   for(; index < n-3; index+=4) {
      __m256d x_index  = _mm256_load_pd(&x[index]);
      __m256d ea_index = _mm256_load_pd(&e->a[index]);
      xea[index] = _mm256_add_pd(x_index, ea_index);
   }
   __m128d xea_rest[3];
   for(; index < n; ++index) {
      __m128d x_index  = _mm_load_sd(&x[index]);
      __m128d ea_index = _mm_load_sd(&e->a[index]);
      xea_rest[index] = _mm_add_sd(x_index, ea_index);
   }

   __m128d cn = _mm128_set_sd(-1.0);
   // Vectorize inner loop in order to profit from contiguous row-major storage
   for(int i=0;i<n;i++) {
      const FT* Ai = Ellipsoid_get_Ai(e,i);
      __m256d ci_pd = _mm256_set_pd(0.0);
      j = 0;
      for(; j<n-3; j+=4) {
        __m256d Aij = _mm256_load_pd(&Ai[j]);
        ci_pd = _mm_fmadd_sd(Aij, xea[j], ci);
      }
      // Reduce and add
      __m128d ci_sd = _mm_set_sd(0.0);
      for(; j<n; j++) {
        __m128d Aij = _mm_load_sd(&Ai[j]);
        ci_sd = _mm_fmadd_sd(Aij, xea[j], ci);
      }
      _mm_store_sd(&c[i], ci_sd);
      cn = _mm_fmadd_sd(ci, xea[i], cn);
   } 
   _mm_store_sd(&c[n], cn);
}*/