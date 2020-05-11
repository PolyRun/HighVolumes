#include "intersectCoord.h"


void Ellipsoid_intersectCoord_cached_reord(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache) {
   Ellipsoid* e = (Ellipsoid*)o;
   const int n = e->n;
   FT* Az_c = (FT*)cache;

   FT a = Ellipsoid_get_a(e,d,d);

   const FT aInv = 0.5/a;

   FT b = 2.0*Az_c[d];
   FT c = Az_c[n];
   
   // find t:
   const FT det = b*b - 4.0*a*c;
   assert(det >= 0);
   const FT sqrtDet = sqrt(det);

   *t0 = (-b - sqrtDet) * aInv;
   *t1 = (-b + sqrtDet) * aInv;
}
