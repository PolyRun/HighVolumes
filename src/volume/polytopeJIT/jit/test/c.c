#include <immintrin.h>
#include <stdbool.h>
#include <stdio.h>

bool inside(const double* x) {
   //return x[1] == 1.0;
   if(x[0]*0.1 + x[1]*1.1 + x[6]*4.2 > 3.76) {return false;}
   if(x[4]*0.5 + x[7]*1.2 + x[8]*5.2 > 3.87) {return false;}
   if(x[3]*0.7 + x[9]*1.4 + x[11]*11.2 > 3.47) {return false;}
   return true;
}


void intersect(const double* x, const double* d, double *t0, double *t1) {
   // ------------- test
   //*t0 = *x;
   //*t1 = *d;
   //return;

   double t00 = -100000.0;
   double t11 = +100000.0;
   
   // for each:
   
   // 2x dot product: d*a, x*a
   
   const double da = 0.1 * d[0] + 9.7*d[1] + 4.4*d[2];
   const double xa = 0.11 * x[0] + 9.71*x[1] + 4.41*x[2];

   // check if d*a is parallel ==0: jump to next - or make sure t is zero?
   double t = (8.5 - xa) / da;
   
   // test:
   __m256d tt = _mm256_set1_pd(t);
   __m256d bb = _mm256_set1_pd(d[2]);
   __m256d cc = _mm256_set1_pd(d[3]);
   __m256d dd = _mm256_max_pd(bb,cc);
   __m256d ee = _mm256_min_pd(dd,tt);
   t = ee[0];

   t00 += t;
   t11 += t;
   // sub, div
   // cond, min, max

   // end
   *t0 = t00;
   *t1 = t11;
}

void intersectCoord_2pack(const double*x,const double*y, double *t0, double*t1) {
   __m128d t00 = _mm_set_pd(-1e10,-1e10);
   
   {
      __m128d xx = _mm_loadu_pd(x+8);
      __m128d yy = _mm_set_pd(0.1,2.0);
      __m128d t = _mm_mul_pd(xx,yy);
      t00 = _mm_max_pd(t00,t);
   }
   {
      __m128d xx = _mm_loadu_pd(x+24);
      __m128d yy = _mm_set_pd(3.3,2.8);
      __m128d t = _mm_mul_pd(xx,yy);
      t00 = _mm_max_pd(t00,t);
   }
   {
      __m128d xx = _mm_loadu_pd(x+40);
      __m128d yy = _mm_set_pd(2.1,7.9);
      __m128d t = _mm_mul_pd(xx,yy);
      t00 = _mm_max_pd(t00,t);
   }
   
   __m128d t00tmp = _mm_permute_pd(t00,0b0101);
   t00 = _mm_max_pd(t00,t00tmp);
   *t0 = t00[0];
}

typedef struct TTT {
   __m256d a;
   __m256d b;
} TTT;

TTT setTTT(double a, double b) {
   TTT ttt;
   ttt.a = _mm256_set1_pd(a);
   ttt.b = _mm256_set1_pd(b);
   return ttt;
}

int main () {
   __m256d a = _mm256_set_pd(4.0, 3.0, 2.0, 1.0);
   __m256d b = _mm256_set_pd(4.0, 3.0, 2.0, 1.1);
   
   __m256d cg = _mm256_cmp_pd(a,b, _CMP_GT_OQ);
   __m256d cl = _mm256_cmp_pd(a,b, _CMP_LT_OQ);
   
   __m256d aa = _mm256_blendv_pd(a,b,cg);
   __m256d bb = _mm256_blendv_pd(a,b,cl);
   
   printf("hello\n");
   printf("hello %f\n",cg[0]);
   printf("hello %f\n",cl[0]);
   printf("hello %f\n",aa[0]);
   printf("hello %f\n",bb[0]);
}
