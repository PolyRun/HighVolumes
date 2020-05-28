#include "linalg.h"

Ball_intersect_f_t Ball_intersect = Ball_intersect_ref;
void Ball_intersect_ref(const int n, const FT r, const FT* x, const FT* d, FT* t0, FT* t1) {
   // y = x + d*t
   // y^2 = r^2
   //
   // x^2 - r^2 + 2x*d*t + d^2*t^2 = 0
   
   FT x2 = squaredNorm(x,n);
   FT d2 = squaredNorm(d,n);
   FT xd = dotProduct(x,d,n);

   FT a = d2;
   FT ainv = 1.0 / a;
   FT b = 2.0*xd;
   FT c = x2 - r*r;

   FT detSqrt = sqrt(b*b - 4.0*a*c);
   
   *t1 = (-b + detSqrt) * 0.5 * ainv;
   *t0 = (-b - detSqrt) * 0.5 * ainv;
}


Ball_intersectCoord_f_t Ball_intersectCoord = Ball_intersectCoord_ref;
void Ball_intersectCoord_ref(const int n, const FT r, const FT* x, const int d, FT* t0, FT* t1) {
   
   FT x2 = squaredNorm(x,n);
   const FT d2 = 1.0;
   FT xd = x[d]; // dot product with unit vector dim d

   const FT a = d2;
   const FT ainv = 1.0 / a;
   FT b = 2.0*xd;
   FT c = x2 - r*r;

   FT detSqrt = sqrt(b*b - 4.0*a*c);
   
   *t1 = (-b + detSqrt) * 0.5 * ainv;
   *t0 = (-b - detSqrt) * 0.5 * ainv;
}

FT Ball_volume(const int n, const FT r) {
   FT rn = pow(r,n);
   
   int nhalf = n / 2;
   bool odd = (n % 2) == 1;
   
   FT pinhalf = pow(M_PI,nhalf);

   FT i = nhalf;
   if(odd) {i+=0.5;}
   
   FT fact = 1.0;
   while(i>0) {fact*=i; i-=1.0;}
   
   return pinhalf * rn / fact;
}

int ceil_cache(const int n, const int b) {
   return (n*b + 31) / 32 * 32;
}


squaredNorm_cached_f_t squaredNorm_cached = squaredNorm_cached_ref;
squaredNorm_cached_reset_f_t squaredNorm_cached_reset = squaredNorm_cached_reset_ref;
squaredNorm_cached_update_f_t squaredNorm_cached_update = squaredNorm_cached_update_ref;


FT squaredNorm_cached_ref(const FT* v, const int n, const FT* cache) {
   // not cached
   return squaredNorm(v,n);
}
void squaredNorm_cached_reset_ref(const FT* v, const int n, FT* cache) {
   // not cached
}
void squaredNorm_cached_update_ref(const FT* v, const int d, const FT dx, const int n, FT* cache) {
   // not cached
}

FT squaredNorm_cached_refc(const FT* v, const int n, const FT* cache) {
   return *cache;
}
void squaredNorm_cached_reset_refc(const FT* v, const int n, FT* cache) {
   *cache = squaredNorm(v,n);
}
void squaredNorm_cached_update_refc(const FT* v, const int d, const FT dx, const int n, FT* cache) {
   FT c = *cache;
   FT x = v[d];
   FT y = x - dx;
   *cache = c - y*y + x*x; // subtract old square, add new square
}

Ball_intersectCoord_cached_f_t Ball_intersectCoord_cached = Ball_intersectCoord_cached_ref;

FTpair Ball_intersectCoord_cached_ref(const int n, const FT r, const FT* x,const int d, FT* cache) {
   FT x2 = squaredNorm_cached(x,n,cache);
   const FT d2 = 1.0;
   FT xd = x[d]; // dot product with unit vector dim d

   const FT a = d2;
   const FT ainv = 1.0 / a;
   FT b = 2.0*xd;
   FT c = x2 - r*r;

   FT detSqrt = sqrt(b*b - 4.0*a*c);
   
   FTpair tp;
   tp.t1 = (-b + detSqrt) * 0.5 * ainv;
   tp.t0 = (-b - detSqrt) * 0.5 * ainv;
   return tp;
}


// ----------------------------------4/8 sets:

__m256d squaredNorm_cached4(const FT* v, const int n, const FT* cache) {
   return _mm256_load_pd(cache);
}
void squaredNorm_cached4_reset(const FT* v, const int n, FT* cache) {
   __m256d acc = _mm256_set1_pd(0.0);
   // TODO: make more efficient!
   for(int i=0;i<n;i++) {
      __m256d vi = _mm256_load_pd(v+4*i);
      acc = _mm256_fmadd_pd(vi,vi,acc);
   }
   _mm256_store_pd(cache,acc);
}
void squaredNorm_cached4_update(const FT* v, const int d, const __m256d dx, const int n, FT* cache) {
   __m256d c = _mm256_load_pd(cache);
   __m256d x = _mm256_load_pd(v+4*d);
   __m256d y = _mm256_sub_pd(x,dx);
   __m256d y2 = _mm256_mul_pd(y,y);
   c = _mm256_fmadd_pd(x,x,c);
   c = _mm256_sub_pd(c,y2);
   _mm256_store_pd(cache,c);
}

FTset8 squaredNorm_cached8(const FT* v, const int n, const FT* cache) {
   __m256d c0 = _mm256_load_pd(cache);
   __m256d c1 = _mm256_load_pd(cache+4);
   FTset8 res = {c0,c1};
   return res;
}
void squaredNorm_cached8_reset(const FT* v, const int n, FT* cache) {
   __m256d acc0 = _mm256_set1_pd(0.0);
   __m256d acc1 = _mm256_set1_pd(0.0);
   // TODO: make more efficient!
   for(int i=0;i<n;i++) {
      __m256d vi0 = _mm256_load_pd(v+8*i);
      __m256d vi1 = _mm256_load_pd(v+8*i+4);
      acc0 = _mm256_fmadd_pd(vi0,vi0,acc0);
      acc1 = _mm256_fmadd_pd(vi1,vi1,acc1);
   }
   _mm256_store_pd(cache,acc0);
   _mm256_store_pd(cache+4,acc1);
}  

void squaredNorm_cached8_update(const FT* v, const int d, const FTset8 dx, const int n, FT* cache) {
   __m256d c0 = _mm256_load_pd(cache);
   __m256d c1 = _mm256_load_pd(cache+4);

   __m256d x0 = _mm256_load_pd(v+8*d);
   __m256d x1 = _mm256_load_pd(v+8*d+4);

   __m256d y0 = _mm256_sub_pd(x0,dx.set0);
   __m256d y1 = _mm256_sub_pd(x1,dx.set1);
   
   __m256d y20 = _mm256_mul_pd(y0,y0);
   __m256d y21 = _mm256_mul_pd(y1,y1);
   
   c0 = _mm256_fmadd_pd(x0,x0,c0);
   c1 = _mm256_fmadd_pd(x1,x1,c1);
   
   c0 = _mm256_sub_pd(c0,y20);
   c1 = _mm256_sub_pd(c1,y21);
   
   _mm256_store_pd(cache,c0);
   _mm256_store_pd(cache+4,c1);
}




