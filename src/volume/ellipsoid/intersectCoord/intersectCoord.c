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
   const FT det = bsquared - a*c;
   assert(det >= 0);
   const FT sqrtDet = sqrt(det);

   *t0 = (-b - sqrtDet) * aInv;
   *t1 = (-b + sqrtDet) * aInv;
   
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
   const FT det = bsquared - a*c; // Issued at 8(12), available at 16(12),    fNmadd here -4
   assert(det >= 0);
   const FT sqrtDet = sqrt(det); //Issued at 16, available at 34

   FT tmp = -b * aInv; // Issued at 14, available at 18

   *t0 = tmp - sqrtDet * aInv; // Issued at 34(38), available at 42(38),     fNmadd here -4
   *t1 = tmp + sqrtDet * aInv; // Issued at 34(38), available at 42(38),     fmadd here  -4
}