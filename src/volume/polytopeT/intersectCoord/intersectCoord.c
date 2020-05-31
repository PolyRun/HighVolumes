#include "intersectCoord.h"
#include <stdio.h>

void PolytopeT_intersectCoord_cached_ref(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache) {
   const PolytopeT* p = (PolytopeT*)o;
   const int n = p->n;
   const int m = p->m;
   FT* Aix = (FT*)cache;
   
   FT t00 = -FT_MAX;// tmp variables for t0, t1
   FT t11 = FT_MAX;

   for(int i=0; i<m; i++) {
      const FT b = PolytopeT_get_b(p, i);
      const FT dai = PolytopeT_get_a(p,i,d); // dot product with unit vector dim d
      
      if(dai <= FT_EPS && -dai <= FT_EPS) {continue;} // orthogonal
      
      const FT aix = Aix[i];
      FT t = (b - aix) / dai;
      
      if(dai < 0.0) {
         t00 = (t00>t)?t00:t; // max
      } else {
         t11 = (t11<t)?t11:t; // min
      }
   }
   
   // return:
   *t0 = t00;
   *t1 = t11;
}

void PolytopeT_intersectCoord_cached_nc1(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache) {
   const PolytopeT* p = (PolytopeT*)o;
   const int n = p->n;
   const int m = p->m;
   FT* Aix = (FT*)cache;
   
   FT t00 = -FT_MAX;// tmp variables for t0, t1
   FT t11 = FT_MAX;

   for(int i=0; i<m; i++) {
      //const FT* ai = PolytopeT_get_Ai(p,i);
      const FT b = PolytopeT_get_b(p, i);
      const FT dai = PolytopeT_get_a(p,i,d); // dot product with unit vector dim d
      
      const FT fwd = (dai > FT_EPS);
      const FT bwd = (dai < -FT_EPS);
      //printf("fwd: %.12f, bwd: %.12f\n",fwd,bwd);
      
      const FT aix = Aix[i];
      FT t = (b - aix) / dai;
      
      t00 = bwd*((t00>t)?t00:t) + (1-bwd)*t00;
      t11 = fwd*((t11<t)?t11:t) + (1-fwd)*t11;
   }
   
   // return:
   *t0 = t00;
   *t1 = t11;
}


void PolytopeT_intersectCoord_cached_b_ref(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache) {
   const PolytopeT* p = (PolytopeT*)o;
   const int n = p->n;
   const int m = p->m;
   FT* cc = (FT*)cache;
   
   FT t00 = -FT_MAX;// tmp variables for t0, t1
   FT t11 = FT_MAX;

   for(int i=0; i<m; i++) {
      //const FT* ai = PolytopeT_get_Ai(p,i);
      //const FT b = PolytopeT_get_b(p, i); MB: don't need this anymore!
      const FT dai = PolytopeT_get_a(p,i,d); // dot product with unit vector dim d
      
      if(dai <= FT_EPS && -dai <= FT_EPS) {continue;} // orthogonal
      
      FT t = cc[i] / dai; // cc[i] = (bi - aix)
      
      if(dai < 0.0) {
         t00 = (t00>t)?t00:t; // max
      } else {
         t11 = (t11<t)?t11:t; // min
      }
   }
   
   // return:
   *t0 = t00;
   *t1 = t11;
}

void PolytopeT_intersectCoord_cached_b_inv_ref(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache) {
   const PolytopeT* p = (PolytopeT*)o;
   const int n = p->n;
   const int m = p->m;
   FT* cc = (FT*)cache;
   
   __m128d zeros = _mm_set_sd(0.0);

   __m128d vt00 = _mm_set_sd(-FT_MAX);
   __m128d vt11 = _mm_set_sd(FT_MAX);

   for(int i=0; i<m; i++) {
      const FT Ainv = PolytopeT_get_aInv(p,i,d); // dot product with unit vector dim d
      
      __m128d vAinv = _mm_set_sd(Ainv);
      __m128d vCache = _mm_set_sd(cc[i]);
      __m128d mask1 = _mm_cmp_sd(vAinv, zeros, _CMP_LT_OS); // towards max
      __m128d mask2 = _mm_cmp_sd(vAinv, zeros, _CMP_GT_OS); // towards min
      __m128d vt = _mm_mul_sd(vAinv,vCache);
      __m128d t00_tmp = _mm_blendv_pd(vt00, vt, mask1); // towards max
      __m128d t11_tmp = _mm_blendv_pd(vt11, vt, mask2); // towards min
      //printf("vAinv %lf %lf %lf \n",vAinv[0],mask1[0],mask2[0]);
      //printf("  %lf %lf %lf\n",vt[0],t00_tmp[0],t11_tmp[0]);
      vt00 = _mm_max_sd(vt00,t00_tmp);
      vt11 = _mm_min_sd(vt11,t11_tmp);

   }
   //printf("res: %lf %lf \n",vt00[0],vt11[0]);
   
   // return:
   *t0 = vt00[0];
   *t1 = vt11[0];
}

void PolytopeT_intersectCoord_cached_b_inv_vec(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache) {
   const PolytopeT* p = (PolytopeT*)o;
   const int n = p->n;
   const int m = p->m;
   FT* cc = (FT*)cache;
   
   //PolytopeT_print(p);

   FT *Ainv_slice = p->Ainv + (p->line * d);
   
   __m256d zeros = _mm256_set1_pd(0.0);

   __m256d vt00 = _mm256_set1_pd(-FT_MAX);
   __m256d vt11 = _mm256_set1_pd(FT_MAX);
   
   for(int i=0; i<m; i+=4) {
      
      __m256d vAinv = _mm256_load_pd(Ainv_slice + i);
      __m256d vCache = _mm256_load_pd(cc + i);
      __m256d mask1 = _mm256_cmp_pd(vAinv, zeros, _CMP_LT_OS); // towards max
      __m256d mask2 = _mm256_cmp_pd(vAinv, zeros, _CMP_GT_OS); // towards min
      __m256d vt = _mm256_mul_pd(vAinv,vCache);
      __m256d t00_tmp = _mm256_blendv_pd(vt00, vt, mask1); // towards max
      __m256d t11_tmp = _mm256_blendv_pd(vt11, vt, mask2); // towards min
      vt00 = _mm256_max_pd(vt00,t00_tmp);
      vt11 = _mm256_min_pd(vt11,t11_tmp);
      //printf("--0  : %lf %lf \n",vAinv[0], vCache[0]);
      //printf("--1  : %lf %lf \n",vAinv[1], vCache[1]);
      //printf("--2  : %lf %lf \n",vAinv[2], vCache[2]);
      //printf("--3  : %lf %lf \n",vAinv[3], vCache[3]);
   }
   //printf("res: %lf %lf \n",vt00[0],vt11[0]);
   
   FT t0_a = vt00[0];
   FT t0_b = (vt00[1] > t0_a) ? vt00[1] : t0_a;
   FT t0_c = (vt00[2] > t0_b) ? vt00[2] : t0_b;
   FT t0_d = (vt00[3] > t0_c) ? vt00[3] : t0_c;

   FT t1_a = vt11[0];
   FT t1_b = (vt11[1] < t1_a) ? vt11[1] : t1_a;
   FT t1_c = (vt11[2] < t1_b) ? vt11[2] : t1_b;
   FT t1_d = (vt11[3] < t1_c) ? vt11[3] : t1_c;
   
   //printf("res: %lf %lf \n",t0_d, t1_d);
   //printf("0  : %lf %lf \n",vt00[0], vt11[0]);
   //printf("1  : %lf %lf \n",vt00[1], vt11[1]);
   //printf("2  : %lf %lf \n",vt00[2], vt11[2]);
   //printf("3  : %lf %lf \n",vt00[3], vt11[3]);

   // return:
   *t0 = t0_d;
   *t1 = t1_d;
}


static inline void inv_loop(const FT *Ainv_slice, const FT *cc, const int i, const __m256d zeros, __m256d *vt00, __m256d *vt11){
          
    __m256d vAinv = _mm256_load_pd(Ainv_slice + i);
    __m256d vCache = _mm256_load_pd(cc + i);
    __m256d mask1 = _mm256_cmp_pd(vAinv, zeros, _CMP_LT_OS); // towards max
    __m256d mask2 = _mm256_cmp_pd(vAinv, zeros, _CMP_GT_OS); // towards min
    __m256d vt = _mm256_mul_pd(vAinv,vCache);
    __m256d t00_tmp = _mm256_blendv_pd(*vt00, vt, mask1); // towards max
    __m256d t11_tmp = _mm256_blendv_pd(*vt11, vt, mask2); // towards min
    *vt00 = _mm256_max_pd(*vt00,t00_tmp);
    *vt11 = _mm256_min_pd(*vt11,t11_tmp);

}


void PolytopeT_intersectCoord_cached_b_inv_vec_inl(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache) {
   const PolytopeT* p = (PolytopeT*)o;
   const int n = p->n;
   const int m = p->m;
   FT* cc = (FT*)cache;
   
   //PolytopeT_print(p);

   FT *Ainv_slice = p->Ainv + (p->line * d);
   
   __m256d zeros = _mm256_set1_pd(0.0);

   __m256d vt001 = _mm256_set1_pd(-FT_MAX);
   __m256d vt111 = _mm256_set1_pd(FT_MAX);
   __m256d vt002 = _mm256_set1_pd(-FT_MAX);
   __m256d vt112 = _mm256_set1_pd(FT_MAX);

   int i = 0;
   for(; i<m-7; i+=8) {
       inv_loop(Ainv_slice, cc, i, zeros, &vt001, &vt111);
       inv_loop(Ainv_slice, cc, i+4, zeros, &vt002, &vt112);
   }
   for (; i < m; i+=4){
       inv_loop(Ainv_slice, cc, i, zeros, &vt001, &vt111);
   }
   __m256d vt00 = _mm256_max_pd(vt001, vt002);
   __m256d vt11 = _mm256_min_pd(vt111, vt112);
   //printf("res: %lf %lf \n",vt00[0],vt11[0]);
   
   FT t0_a = vt00[0];
   FT t0_b = (vt00[1] > t0_a) ? vt00[1] : t0_a;
   FT t0_c = (vt00[2] > t0_b) ? vt00[2] : t0_b;
   *t0 = (vt00[3] > t0_c) ? vt00[3] : t0_c;

   FT t1_a = vt11[0];
   FT t1_b = (vt11[1] < t1_a) ? vt11[1] : t1_a;
   FT t1_c = (vt11[2] < t1_b) ? vt11[2] : t1_b;
   *t1 = (vt11[3] < t1_c) ? vt11[3] : t1_c;
   
}



void PolytopeT_intersectCoord_cached_b_vec(const void *p, const FT *x, const int d,
                                           FT *t0_out, FT *t1_out, void *cache) {

   const PolytopeT *poly = (PolytopeT*) p;
   const int dims = poly->n;
   const int constraints = poly->m;

   FT *b_sub_Aix_slice = (FT*) cache;
   FT *Aid_slice = poly->A + (poly->line * d);

   __m256d t0_vec_0 = _mm256_set1_pd(-FT_MAX);
   __m256d t1_vec_0 = _mm256_set1_pd(FT_MAX);

   __m256d ft_neg_vec = _mm256_set1_pd(-FT_EPS);
   __m256d ft_pos_vec = _mm256_set1_pd(FT_EPS);

   int i = 0;
   for (; i < constraints - 3; i += 4) {

      __m256d b_sub_Aix_vec_0 = _mm256_load_pd(b_sub_Aix_slice + i);
      __m256d Aid_vec_0 = _mm256_load_pd(Aid_slice + i);

      __m256d t_vec_0 = _mm256_div_pd(b_sub_Aix_vec_0, Aid_vec_0);

      __m256d neg_mask = _mm256_cmp_pd(Aid_vec_0, ft_neg_vec, _CMP_LE_OQ);
      __m256d pos_mask = _mm256_cmp_pd(Aid_vec_0, ft_pos_vec, _CMP_GE_OQ);

      __m256d tmp_t0_vec_0 = _mm256_blendv_pd(t0_vec_0, t_vec_0, neg_mask);
      __m256d tmp_t1_vec_0 = _mm256_blendv_pd(t1_vec_0, t_vec_0, pos_mask);

      t0_vec_0 = _mm256_max_pd(t0_vec_0, tmp_t0_vec_0);
      t1_vec_0 = _mm256_min_pd(t1_vec_0, tmp_t1_vec_0);
   }

   FT t0;
   FT t1;

   t0 = t0_vec_0[0];
   t0 = (t0_vec_0[1] > t0) ? t0_vec_0[1] : t0;
   t0 = (t0_vec_0[2] > t0) ? t0_vec_0[2] : t0;
   t0 = (t0_vec_0[3] > t0) ? t0_vec_0[3] : t0;

   t1 = t1_vec_0[0];
   t1 = (t1_vec_0[1] < t1) ? t1_vec_0[1] : t1;
   t1 = (t1_vec_0[2] < t1) ? t1_vec_0[2] : t1;
   t1 = (t1_vec_0[3] < t1) ? t1_vec_0[3] : t1;

   for (; i < constraints; i++) {

      FT b_sub_Aix = b_sub_Aix_slice[i];
      FT Aid = Aid_slice[i];

      FT t = b_sub_Aix / Aid;

      if (Aid < 0.0) {
         t0 = (t0 > t) ? t0 : t;
      } else {
         t1 = (t1 < t) ? t1 : t;
      }
   }

   *t0_out = t0;
   *t1_out = t1;
}


/*
 NOTE: this implementation assumes that the invalid elements of A are set to 0
 as it doesn't outperform cached_b_vec it can just be deleted if this assumption should pose a problem at some point...
 */
void PolytopeT_intersectCoord_cached_b_vec2(const void *p, const FT *x, const int d,
                                           FT *t0, FT *t1, void *cache) {

   const PolytopeT *poly = (PolytopeT*) p;
   const int n = poly->n;
   // now m is divisible by 4
   // moreover, we don't care about divisions by 0 as they give +-inf which doesn't affect our max/min computation
   // also invalid elements should be set to 0 by initialization and are then never touched again, right?
   const int m = poly->line / sizeof(FT);

   FT *b_Aix = (FT*) cache;
   FT *Aid = poly->A + (poly->line * d);

   __m256d t0_vec_0 = _mm256_set1_pd(-FT_MAX);
   __m256d t1_vec_0 = _mm256_set1_pd(FT_MAX);
   __m256d zeros = _mm256_set1_pd(0.0);

   for (int i = 0; i < m - 3; i += 4) {

      __m256d b_sub_Aix_v = _mm256_load_pd(b_Aix + i);
      __m256d Aid_v = _mm256_load_pd(Aid + i);

      __m256d t = _mm256_div_pd(b_sub_Aix_v, Aid_v);

      __m256d n_mask = _mm256_cmp_pd(t, zeros, _CMP_LT_OS);
      __m256d p_mask = _mm256_cmp_pd(t, zeros, _CMP_GE_OS);

      __m256d tmp_t0 = _mm256_blendv_pd(t0_vec_0, t, n_mask);
      __m256d tmp_t1 = _mm256_blendv_pd(t1_vec_0, t, p_mask);

      t0_vec_0 = _mm256_max_pd(t0_vec_0, tmp_t0);
      t1_vec_0 = _mm256_min_pd(t1_vec_0, tmp_t1);
   }
   
   FT t0t0t0 = (t0_vec_0[1] > t0_vec_0[0]) ? t0_vec_0[1] : t0_vec_0[0];
   FT t0t0 = (t0_vec_0[2] > t0t0t0) ? t0_vec_0[2] : t0t0t0;
   *t0 = (t0_vec_0[3] > t0t0) ? t0_vec_0[3] : t0t0;

   FT t1t1t1 = (t1_vec_0[1] < t1_vec_0[0]) ? t1_vec_0[1] : t1_vec_0[0];
   FT t1t1 = (t1_vec_0[2] < t1t1t1) ? t1_vec_0[2] : t1t1t1;
   *t1 = (t1_vec_0[3] < t1t1) ? t1_vec_0[3] : t1t1;
}


static inline void loop(const FT *b_Aix, const FT *Aid, const int i, const __m256d zeros, __m256d *t0_vec, __m256d *t1_vec){
    __m256d b_sub_Aix_v = _mm256_load_pd(b_Aix + i);
    __m256d Aid_v = _mm256_load_pd(Aid + i);

    __m256d t = _mm256_div_pd(b_sub_Aix_v, Aid_v);

    __m256d n_mask = _mm256_cmp_pd(t, zeros, _CMP_LT_OS);
    __m256d p_mask = _mm256_cmp_pd(t, zeros, _CMP_GE_OS);

    __m256d tmp_t0 = _mm256_blendv_pd(*t0_vec, t, n_mask);
    __m256d tmp_t1 = _mm256_blendv_pd(*t1_vec, t, p_mask);

    *t0_vec = _mm256_max_pd(*t0_vec, tmp_t0);
    *t1_vec = _mm256_min_pd(*t1_vec, tmp_t1);

}


void PolytopeT_intersectCoord_cached_b_vec_inl(const void *p, const FT *x, const int d,
                                           FT *t0, FT *t1, void *cache) {

   const PolytopeT *poly = (PolytopeT*) p;
   const int n = poly->n;
   // now m is divisible by 4
   // moreover, we don't care about divisions by 0 as they give +-inf which doesn't affect our max/min computation
   // also invalid elements should be set to 0 by initialization and are then never touched again, right?
   const int m = poly->line / sizeof(FT);

   FT *b_Aix = (FT*) cache;
   FT *Aid = poly->A + (poly->line * d);

   __m256d t0_vec_0 = _mm256_set1_pd(-FT_MAX);
   __m256d t1_vec_0 = _mm256_set1_pd(FT_MAX);
   __m256d t0_vec_1 = _mm256_set1_pd(-FT_MAX);
   __m256d t1_vec_1 = _mm256_set1_pd(FT_MAX);
   __m256d t0_vec_2 = _mm256_set1_pd(-FT_MAX);
   __m256d t1_vec_2 = _mm256_set1_pd(FT_MAX);
   __m256d t0_vec_3 = _mm256_set1_pd(-FT_MAX);
   __m256d t1_vec_3 = _mm256_set1_pd(FT_MAX);
   __m256d zeros = _mm256_set1_pd(0.0);

   int i = 0;
   for (; i < m - 15; i += 16) {

       loop(b_Aix, Aid, i, zeros, &t0_vec_0, &t1_vec_0);
       loop(b_Aix, Aid, i+4, zeros, &t0_vec_1, &t1_vec_1);
       loop(b_Aix, Aid, i+8, zeros, &t0_vec_2, &t1_vec_2);
       loop(b_Aix, Aid, i+12, zeros, &t0_vec_3, &t1_vec_3);
       
   }
   for (; i < m-3; i+=4){
       loop(b_Aix, Aid, i, zeros, &t0_vec_0, &t1_vec_0);
   }

   
   __m256d t00 = _mm256_max_pd(t0_vec_0, t0_vec_1);
   __m256d t01 = _mm256_max_pd(t0_vec_2, t0_vec_3);
   __m256d t000 = _mm256_max_pd(t00, t01);
   __m256d t10 = _mm256_min_pd(t1_vec_0, t1_vec_1);
   __m256d t11 = _mm256_min_pd(t1_vec_2, t1_vec_3);
   __m256d t111 = _mm256_min_pd(t10, t11);
   
   
   FT t0t0t0 = (t000[1] > t000[0]) ? t000[1] : t000[0];
   FT t0t0 = (t000[2] > t0t0t0) ? t000[2] : t0t0t0;
   *t0 = (t000[3] > t0t0) ? t000[3] : t0t0;

   FT t1t1t1 = (t111[1] < t111[0]) ? t111[1] : t111[0];
   FT t1t1 = (t111[2] < t1t1t1) ? t111[2] : t1t1t1;
   *t1 = (t111[3] < t1t1) ? t111[3] : t1t1;
}

void PolytopeT_cacheReset_b_ref(const void* o, const FT* x, void* cache) {
   const PolytopeT* p = (PolytopeT*)o;
   FT* c = (FT*)cache; // set c[i] = bi - Ai*x
   const int n = p->n;
   const int m = p->m;
   for(int i=0; i<m; i++) {
      FT dot = PolytopeT_get_b(p,i); // watch the b!
      for(int j=0;j<n;j++) {
         dot -= x[j] * PolytopeT_get_a(p,i,j);  // watch the minus!
      }
      c[i] = dot;
   }
}

void PolytopeT_cacheReset_b_vec(const void *p, const FT *x, void *cache) {

   const PolytopeT *poly = (PolytopeT*) p;
   FT *cache_slice = (FT*) cache;
   FT *A_slice =poly->A;
   FT *b_slice = poly->b;

   const int dims = poly->n;
   const int constraints = poly->m;

   // Initializing the cache with b[i], so we can later subtract Ai dot x from it
   int i;
   for (i = 0; i < constraints - 4; i += 4) {
      __m256d b_vec = _mm256_load_pd(b_slice + i);
      _mm256_store_pd(cache_slice + i, b_vec);
   }
   for (; i < constraints; i++) {
      cache_slice[i] = b_slice[i];
   }

   // Accessing row-major (i.e. memory-aligned since the polytope is transposed)
   int j;
   for (j = 0; j < dims; j++) {

      FT *Aij_slice = A_slice + (poly->line * j);
      FT x_j = x[j];
      __m256d x_j_vec = _mm256_set1_pd(x_j);

      for (i = 0; i < constraints - 4; i += 4) {
         __m256d cache_vec = _mm256_load_pd(cache_slice + i);
         __m256d A_vec = _mm256_load_pd(Aij_slice + i);

         __m256d A_dot_x_vec = _mm256_mul_pd(A_vec, x_j_vec);
         cache_vec = _mm256_sub_pd(cache_vec, A_dot_x_vec);

         _mm256_store_pd(cache_slice + i, cache_vec);
      }

      for (; i < constraints; i++) {
         cache_slice[i] -= x_j * Aij_slice[i];
      }

   }
}

void PolytopeT_cacheUpdateCoord_b_ref(const void* o, const int d, const FT dx, void* cache) {
   const PolytopeT* p = (PolytopeT*)o;
   const int m = p->m;
   FT* c = (FT*)cache;
   for(int i=0; i<m; i++) {
      c[i] -= dx * PolytopeT_get_a(p,i,d); // watch the minus !
   } 
}

void PolytopeT_cacheUpdateCoord_b_vec(const void *p, const int d, const FT dx, void *cache) {

   const PolytopeT *poly = (PolytopeT*) p;
   const int constraints = poly->m;

   FT *cache_slice = (FT*) cache;
   FT *Aid_slice = poly->A + (poly->line * d);
   
   __m256d dx_vec = _mm256_set1_pd(dx);

   int i;
   for (i = 0; i < constraints; i += 4) {
      
      __m256d Aid_vec = _mm256_load_pd(Aid_slice + i);
      __m256d cache_vec = _mm256_load_pd(cache_slice + i);

      __m256d Aid_mul_dx_vec = _mm256_mul_pd(Aid_vec, dx_vec);
      __m256d update_vec = _mm256_sub_pd(cache_vec, Aid_mul_dx_vec);

      _mm256_store_pd(cache_slice + i, update_vec);
   }

   for (; i < constraints; i++) {
      cache_slice[i] -= dx * Aid_slice[i];
   }
}



void PolytopeT_intersectCoord_vectorized(const void *p, const FT* x,
                                         const int d, FT *t0_out, FT *t1_out, void *cache) {

   const PolytopeT *poly = (PolytopeT*) p;
   const int dims = poly->n;
   const int constraints = poly->m;

   // We will iteratively process slices from b, starting at &b[0]
   FT *b_slice = poly->b;

   // Similarly for the dot product A_i * x
   FT *Aix_slice = (FT*) cache;

   // And the dot product A_i * e_d, where e_d is the d-th unit vector
   FT *Aid_slice = poly->A + (poly->line * d);

   // The same t0 and t1 as in the other implementations
   // i.e. t0 is the lower bound and t1 the upper bound of our segment
   __m256d t0_vec_0, t0_vec_1;
   __m256d t1_vec_0, t1_vec_1;
   t0_vec_0 = _mm256_set1_pd(-FT_MAX);
   t0_vec_1 = _mm256_set1_pd(-FT_MAX);
   t1_vec_0 = _mm256_set1_pd(FT_MAX);
   t1_vec_1 = _mm256_set1_pd(FT_MAX);

   // Constants for FT_EPS and -FT_EPS
   __m256d ft_eps_pos = _mm256_set1_pd(FT_EPS);
   __m256d ft_eps_neg = _mm256_set1_pd(-FT_EPS);

   // Before continuing, we should discuss one design decision.
   // As already seen in the reference implementation, we need to test
   // whether Aid is close to zero, because if so, we don't do anything.
   // Furthermore we need check, whether our new t is negative or positive
   // and update t0 or t1 resp., if t is bigger than t0 or smaller than t1.
   // In total this gives 4 check: 2 for Aid (1 '<=' and 1 '>=') and 

   int i;
   for (i = 0; i < constraints - 8; i += 8) {

      __m256d b_vec_0 = _mm256_load_pd(b_slice + i);
      __m256d b_vec_1 = _mm256_load_pd(b_slice + i + 4);

      __m256d Aix_vec_0 = _mm256_load_pd(Aix_slice + i);
      __m256d Aix_vec_1 = _mm256_load_pd(Aix_slice + i + 4);

      __m256d Aid_vec_0 = _mm256_load_pd(Aid_slice + i);
      __m256d Aid_vec_1 = _mm256_load_pd(Aid_slice + i + 4);

      Aix_vec_0 = _mm256_sub_pd(b_vec_0, Aix_vec_0);
      Aix_vec_1 = _mm256_sub_pd(b_vec_1, Aix_vec_1);

      __m256d t_vec_0 = _mm256_div_pd(Aix_vec_0, Aid_vec_0);
      __m256d t_vec_1 = _mm256_div_pd(Aix_vec_1, Aid_vec_1);

      // t_vec_x now contains what was 't' in the usual algorithms

      __m256d mask_leq_0 = _mm256_cmp_pd(Aid_vec_0, ft_eps_neg, _CMP_LE_OQ);
      __m256d mask_leq_1 = _mm256_cmp_pd(Aid_vec_1, ft_eps_neg, _CMP_LE_OQ);
      
      __m256d mask_geq_0 = _mm256_cmp_pd(Aid_vec_0, ft_eps_pos, _CMP_GE_OQ);
      __m256d mask_geq_1 = _mm256_cmp_pd(Aid_vec_1, ft_eps_pos, _CMP_GE_OQ);

      __m256d nonorthogonal_0 = _mm256_or_pd(mask_leq_0, mask_geq_0);
      __m256d nonorthogonal_1 = _mm256_or_pd(mask_leq_1, mask_geq_1);

      // Now we update t0

      __m256d negative_mask_0 = _mm256_cmp_pd(Aid_vec_0, _mm256_setzero_pd(), _CMP_LE_OS);
      __m256d negative_mask_1 = _mm256_cmp_pd(Aid_vec_1, _mm256_setzero_pd(), _CMP_LE_OS);

      __m256d bigger_0  = _mm256_cmp_pd(t_vec_0, t0_vec_0, _CMP_GE_OQ);
      __m256d bigger_1  = _mm256_cmp_pd(t_vec_1, t0_vec_1, _CMP_GE_OQ);

      __m256d update_mask_0 = _mm256_and_pd(bigger_0, nonorthogonal_0);
      __m256d update_mask_1 = _mm256_and_pd(bigger_1, nonorthogonal_1);

      update_mask_0 = _mm256_and_pd(update_mask_0, negative_mask_0);
      update_mask_1 = _mm256_and_pd(update_mask_1, negative_mask_1);

      t0_vec_0 = _mm256_blendv_pd(t0_vec_0, t_vec_0, update_mask_0);
      t0_vec_1 = _mm256_blendv_pd(t0_vec_1, t_vec_1, update_mask_1);

      // Now we update t1

      __m256d positive_mask_0 = _mm256_cmp_pd(Aid_vec_0, _mm256_setzero_pd(), _CMP_GE_OQ);
      __m256d positive_mask_1 = _mm256_cmp_pd(Aid_vec_1, _mm256_setzero_pd(), _CMP_GE_OQ);

      __m256d smaller_0 = _mm256_cmp_pd(t_vec_0, t1_vec_0, _CMP_LE_OS);
      __m256d smaller_1 = _mm256_cmp_pd(t_vec_1, t1_vec_1, _CMP_LE_OS);

      update_mask_0 = _mm256_and_pd(smaller_0, nonorthogonal_0);
      update_mask_1 = _mm256_and_pd(smaller_1, nonorthogonal_1);

      update_mask_0 = _mm256_and_pd(update_mask_0, positive_mask_0);
      update_mask_1 = _mm256_and_pd(update_mask_1, positive_mask_1);

      t1_vec_0 = _mm256_blendv_pd(t1_vec_0, t_vec_0, update_mask_0);
      t1_vec_1 = _mm256_blendv_pd(t1_vec_1, t_vec_1, update_mask_1);

   }

   // Before we finish the remaining iterations we can combine our t0's and t1's

   t0_vec_0 = _mm256_max_pd(t0_vec_0, t0_vec_1);
   t1_vec_0 = _mm256_min_pd(t1_vec_0, t1_vec_1);

   // Now we finish the remaining iterations that the stepsize of 8 left out

   for (; i < constraints - 4; i += 4) {
      __m256d b_vec   = _mm256_load_pd(b_slice + i);
      __m256d Aix_vec = _mm256_load_pd(Aix_slice + i);
      __m256d Aid_vec = _mm256_load_pd(Aid_slice + i);

      Aix_vec = _mm256_sub_pd(b_vec, Aix_vec);
      __m256d t_vec = _mm256_div_pd(Aix_vec, Aid_vec);

      __m256d mask_leq = _mm256_cmp_pd(Aid_vec, ft_eps_neg, _CMP_LE_OS);
      __m256d mask_geq = _mm256_cmp_pd(Aid_vec, ft_eps_pos, _CMP_GE_OS);

      __m256d nonorthogonal =_mm256_or_pd(mask_leq, mask_geq);

      // Update t0

      __m256d bigger = _mm256_cmp_pd(t_vec, t0_vec_0, _CMP_GE_OS);
      __m256d negative = _mm256_cmp_pd(Aid_vec, _mm256_setzero_pd(), _CMP_LE_OS);
      __m256d update_mask = _mm256_and_pd(bigger, nonorthogonal);
      update_mask = _mm256_and_pd(update_mask, negative);
      t0_vec_0 = _mm256_blendv_pd(t0_vec_0, t_vec, update_mask);

      // Update t1

      __m256d smaller = _mm256_cmp_pd(t_vec, t1_vec_0, _CMP_LE_OS);
      __m256d positive = _mm256_cmp_pd(Aid_vec, _mm256_setzero_pd(), _CMP_GE_OQ);
      update_mask = _mm256_and_pd(smaller, nonorthogonal);
      update_mask = _mm256_and_pd(update_mask, positive);
      t1_vec_0 = _mm256_blendv_pd(t1_vec_0, t_vec, update_mask);

   }

   // Before we finish the last last iterations we conbine the t0's and t1's

   FT t0;
   FT t1;

   t0 = t0_vec_0[0];
   t0 = (t0_vec_0[1] > t0) ? t0_vec_0[1] : t0;
   t0 = (t0_vec_0[2] > t0) ? t0_vec_0[2] : t0;
   t0 = (t0_vec_0[3] > t0) ? t0_vec_0[3] : t0;

   t1 = t1_vec_0[0];
   t1 = (t1_vec_0[1] < t1) ? t1_vec_0[1] : t1;
   t1 = (t1_vec_0[2] < t1) ? t1_vec_0[2] : t1;
   t1 = (t1_vec_0[3] < t1) ? t1_vec_0[3] : t1;

   for (; i < constraints; i++) {

      const FT Aid = Aid_slice[i];

      if (-FT_EPS <= Aid && Aid <= FT_EPS) {
         continue;
      }

      const FT b   = b_slice[i];
      const FT Aix = Aix_slice[i];

      const FT t = (b - Aix) / Aid;

      if (Aid < 0) {
         t0 = (t0 > t) ? t0 : t;
      } else {
         t1 = (t1 < t) ? t1 : t;
      }
   }
   
   *t0_out = t0;
   *t1_out = t1;
}



FTpair4 PolytopeT_intersectCoord4_ref(const void* o, const FT* x, const int d, void* cache) {
   const PolytopeT* p = (PolytopeT*)o;
   const int n = p->n;
   const int m = p->m;
   FT* cc = (FT*)cache;
   
   __m256d zeros = _mm256_set1_pd(0.0);

   __m256d vt00 = _mm256_set1_pd(-FT_MAX);
   __m256d vt11 = _mm256_set1_pd(FT_MAX);

   for(int i=0; i<m; i++) {
      // const FT Ainv = PolytopeT_get_aInv(p,i,d); // dot product with unit vector dim d
      // __m128d vAinv = _mm_set_sd(Ainv);
      __m256d vAinv = _mm256_broadcast_sd(p->Ainv + p->line*d + i);
      //// printf("vAinv[%d]: %lf %lf %lf %lf\n",i,vAinv[0],vAinv[1],vAinv[2],vAinv[3]);

      // __m128d vCache = _mm_set_sd(cc[i]);
      __m256d vCache = _mm256_load_pd(cc+4*i);
      //// printf("vC[%d]: %lf %lf %lf %lf\n",i,vCache[0],vCache[1],vCache[2],vCache[3]);
      // __m128d mask1 = _mm_cmp_sd(vAinv, zeros, _CMP_LT_OS); // towards max
      // __m128d mask2 = _mm_cmp_sd(vAinv, zeros, _CMP_GT_OS); // towards min
      __m256d mask1 = _mm256_cmp_pd(vAinv, zeros, _CMP_LT_OS);
      __m256d mask2 = _mm256_cmp_pd(vAinv, zeros, _CMP_GT_OS);
      // __m128d vt = _mm_mul_sd(vAinv,vCache);
      __m256d vt = _mm256_mul_pd(vAinv,vCache);
      //// printf("vt[%d]: %lf %lf %lf %lf\n",i,vt[0],vt[1],vt[2],vt[3]);
      // __m128d t00_tmp = _mm_blendv_pd(vt00, vt, mask1); // towards max
      // __m128d t11_tmp = _mm_blendv_pd(vt11, vt, mask2); // towards min
      __m256d t00_tmp = _mm256_blendv_pd(vt00, vt, mask1);
      __m256d t11_tmp = _mm256_blendv_pd(vt11, vt, mask2);
      // vt00 = _mm_max_sd(vt00,t00_tmp);
      // vt11 = _mm_min_sd(vt11,t11_tmp);
      vt00 = _mm256_max_pd(vt00,t00_tmp);
      vt11 = _mm256_min_pd(vt11,t11_tmp);
      //// printf("vt00[%d]: %lf %lf %lf %lf\n",i,vt00[0],vt00[1],vt00[2],vt00[3]);
      //// printf("vt11[%d]: %lf %lf %lf %lf\n",i,vt11[0],vt11[1],vt11[2],vt11[3]);
   }
   ////printf("res: %lf %lf \n",vt00[0],vt11[0]);
   //
   //// return:
   //*t0 = vt00[0];
   //*t1 = vt11[0];
   FTpair4 tp;
   tp.low0 = vt00;
   tp.hi0  = vt11;
   //// printf("tp.lo: %lf %lf %lf %lf\n",tp.low0[0],tp.low0[1],tp.low0[2],tp.low0[3]);
   //// printf("tp.hi: %lf %lf %lf %lf\n",tp.hi0[0],tp.hi0[1],tp.hi0[2],tp.hi0[3]);
   return tp;
}
FTpair8 PolytopeT_intersectCoord8_ref(const void* o, const FT* x, const int d, void* cache) {
   const PolytopeT* p = (PolytopeT*)o;
   const int n = p->n;
   const int m = p->m;
   FT* cc = (FT*)cache;
   
   __m256d zeros = _mm256_set1_pd(0.0);

   __m256d vt00 = _mm256_set1_pd(-FT_MAX);
   __m256d vt01 = _mm256_set1_pd(-FT_MAX);
   __m256d vt10 = _mm256_set1_pd(FT_MAX);
   __m256d vt11 = _mm256_set1_pd(FT_MAX);

   for(int i=0; i<m; i++) {
      __m256d vAinv = _mm256_broadcast_sd(p->Ainv + p->line*d + i);
      
      __m256d vCache0 = _mm256_load_pd(cc+8*i);
      __m256d vCache1 = _mm256_load_pd(cc+8*i+4);
      
      __m256d mask1 = _mm256_cmp_pd(vAinv, zeros, _CMP_LT_OS);
      __m256d mask2 = _mm256_cmp_pd(vAinv, zeros, _CMP_GT_OS);
      
      __m256d vt0 = _mm256_mul_pd(vAinv,vCache0);
      __m256d vt1 = _mm256_mul_pd(vAinv,vCache1);
      
      __m256d t00_tmp = _mm256_blendv_pd(vt00, vt0, mask1);
      __m256d t01_tmp = _mm256_blendv_pd(vt01, vt1, mask1);
      __m256d t10_tmp = _mm256_blendv_pd(vt10, vt0, mask2);
      __m256d t11_tmp = _mm256_blendv_pd(vt11, vt1, mask2);
      
      vt00 = _mm256_max_pd(vt00,t00_tmp);
      vt01 = _mm256_max_pd(vt01,t01_tmp);
      vt10 = _mm256_min_pd(vt10,t10_tmp);
      vt11 = _mm256_min_pd(vt11,t11_tmp);
   }
   FTpair8 tp;
   tp.low0 = vt00;
   tp.low1 = vt01;
   tp.hi0  = vt10;
   tp.hi1  = vt11;
   return tp;

}
void PolytopeT_cacheReset4_ref(const void* o, const FT* x, void* cache) {
   const PolytopeT* p = (PolytopeT*)o;
   FT* c = (FT*)cache; // set c[i] = bi - Ai*x
   const int n = p->n;
   const int m = p->m;
   const FT* b = p->b;
   for(int i=0; i<m; i++) {
      // FT dot = PolytopeT_get_b(p,i); // watch the b!
      __m256d dot = _mm256_broadcast_sd(b+i);
      for(int j=0;j<n;j++) {
         // dot -= x[j] * PolytopeT_get_a(p,i,j);  // watch the minus!
         __m256d xj = _mm256_load_pd(x+4*j);
	 __m256d aij = _mm256_broadcast_sd(p->A + p->line*j + i);
	 __m256d xa = _mm256_mul_pd(xj,aij); // can this be optimized?
	 dot = _mm256_sub_pd(dot,xa);
      }
      // c[i] = dot;
      _mm256_store_pd(c+4*i, dot);
      //// printf("c[%d]: %lf %lf %lf %lf\n",i,dot[0],dot[1],dot[2],dot[3]);
   }
}
void PolytopeT_cacheReset8_ref(const void *o, const FT *x, void *cache) {
   const PolytopeT* p = (PolytopeT*)o;
   FT* c = (FT*)cache; // set c[i] = bi - Ai*x
   const int n = p->n;
   const int m = p->m;
   const FT* b = p->b;
   for(int i=0; i<m; i++) {
      __m256d dot0 = _mm256_broadcast_sd(b+i);
      __m256d dot1 = dot0;//_mm256_broadcast_sd(b+i);
      for(int j=0;j<n;j++) {
         __m256d xj0 = _mm256_load_pd(x+8*j);
         __m256d xj1 = _mm256_load_pd(x+8*j+4);
	 __m256d aij = _mm256_broadcast_sd(p->A + p->line*j + i);
	 __m256d xa0 = _mm256_mul_pd(xj0,aij); // can this be optimized?
	 __m256d xa1 = _mm256_mul_pd(xj1,aij); // can this be optimized?
	 dot0 = _mm256_sub_pd(dot0,xa0);
	 dot1 = _mm256_sub_pd(dot1,xa1);
      }
      _mm256_store_pd(c+8*i,   dot0);
      _mm256_store_pd(c+8*i+4, dot1);
   }
}
// void PolytopeT_cacheUpdateCoord4_ref(const void* o, const int d, const __m256d dx, void* cache) {
//    const PolytopeT* p = (PolytopeT*)o;
//    const int m = p->m;
//    FT* c = (FT*)cache;
//    FT* offset = p->A + p->line*d;
//    for(int i=0; i<m; i++) {
//       //FT* cc0 = c+4*i;
//       //__m256d ci = _mm256_load_pd(cc0);
//       __m256d ci = _mm256_load_pd(c);
//       __m256d aid = _mm256_broadcast_sd(offset + i);
//       __m256d tmp = _mm256_mul_pd(dx,aid);
//       __m256d tmp_ = _mm256_sub_pd(ci,tmp);
//       //_mm256_store_pd(cc0, tmp_);
//       _mm256_store_pd(c, tmp_);
//       c +=4;
//    } 
// }
void PolytopeT_cacheUpdateCoord4_ref(const void* o, const int d, const __m256d dx, void* cache) {
   const PolytopeT* p = (PolytopeT*)o;
   const int m = p->m;
   FT* c = (FT*)cache;
   FT* offset = p->A + p->line*d;
   int i=0;
   for(; i<m-3; i+=4) {
      __m256d ci0 = _mm256_load_pd(c);
      __m256d ci1 = _mm256_load_pd(c+4);
      __m256d ci2 = _mm256_load_pd(c+8);
      __m256d ci3 = _mm256_load_pd(c+12);
      __m256d aid0 = _mm256_broadcast_sd(offset + i);
      __m256d aid1 = _mm256_broadcast_sd(offset + i+1);
      __m256d aid2 = _mm256_broadcast_sd(offset + i+2);
      __m256d aid3 = _mm256_broadcast_sd(offset + i+3);
      __m256d tmp0 = _mm256_mul_pd(dx,aid0);
      __m256d tmp1 = _mm256_mul_pd(dx,aid1);
      __m256d tmp2 = _mm256_mul_pd(dx,aid2);
      __m256d tmp3 = _mm256_mul_pd(dx,aid3);
      __m256d tmp_0 = _mm256_sub_pd(ci0,tmp0);
      __m256d tmp_1 = _mm256_sub_pd(ci1,tmp1);
      __m256d tmp_2 = _mm256_sub_pd(ci2,tmp2);
      __m256d tmp_3 = _mm256_sub_pd(ci3,tmp3);
      _mm256_store_pd(c, tmp_0);
      _mm256_store_pd(c+4, tmp_1);
      _mm256_store_pd(c+8, tmp_2);
      _mm256_store_pd(c+12, tmp_3);
      c+=16;
   } 
   for(; i<m; i++) {
      //FT* cc0 = c+4*i;
      //__m256d ci = _mm256_load_pd(cc0);
      __m256d ci = _mm256_load_pd(c);
      __m256d aid = _mm256_broadcast_sd(offset + i);
      __m256d tmp = _mm256_mul_pd(dx,aid);
      __m256d tmp_ = _mm256_sub_pd(ci,tmp);
      //_mm256_store_pd(cc0, tmp_);
      _mm256_store_pd(c, tmp_);
      c +=4;
   } 
}
void PolytopeT_cacheUpdateCoord8_ref(const void* o, const int d, const FTset8 dx, void* cache) {
   const PolytopeT* p = (PolytopeT*)o;
   const int m = p->m;
   FT* c = (FT*)cache;
   FT* offset = p->A + p->line*d;
   __m256d dx0 = dx.set0;
   __m256d dx1 = dx.set1;
   for(int i=0; i<m; i++) {
      FT* cc0 = c+8*i;
      FT* cc1 = cc0+4;
      __m256d ci0 = _mm256_load_pd(cc0);
      __m256d ci1 = _mm256_load_pd(cc1);
      __m256d aid = _mm256_broadcast_sd(offset + i);
      //__m256d tmp0 = _mm256_mul_pd(dx.set0,aid);
      //__m256d tmp1 = _mm256_mul_pd(dx.set1,aid);
      __m256d tmp0 = _mm256_mul_pd(dx0,aid);
      __m256d tmp1 = _mm256_mul_pd(dx1,aid);
      __m256d tmp_0 = _mm256_sub_pd(ci0,tmp0);
      __m256d tmp_1 = _mm256_sub_pd(ci1,tmp1);
      _mm256_store_pd(cc0, tmp_0);
      _mm256_store_pd(cc1, tmp_1);
   } 
}


