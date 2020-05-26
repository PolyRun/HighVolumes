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


