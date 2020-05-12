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

   printf("%f and %f \n", t0, t1);
   
   *t0_out = t0;
   *t1_out = t1;
}
