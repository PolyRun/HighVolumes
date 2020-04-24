#include "intersectCoord.h"

void PolytopeT_intersectCoord_cached_ref(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache) {
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

