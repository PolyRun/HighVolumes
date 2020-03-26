#include "volume.h"


FT dotProduct(const FT* u, const FT* v, const int n) {
   FT sum = 0.0;
   for(int i=0; i<n; i++) {sum+= u[i]*v[i];}
   return sum;
}

Polytope* Polytope_new(int n, int m) {
   Polytope* o = (Polytope*) malloc(sizeof(Polytope));
   o->data = (FT*) malloc(sizeof(FT)*(n+1)*m);
   o->n = n;
   o->m = m;

   return o;
}


void Polytope_free(Polytope* p) {
   free(p->data);
   free(p);
}

FT* Polytope_get_aV(const Polytope* p, int i) {
   return &(p->data[i * (p->n+1)]);
}

bool Polytope_inside(const Polytope* p, const FT* v) {
   for(int i=0; i<p->n; i++) {
      FT sum = 0;
      for(int x=0; x<p->n; x++) { sum+= v[x] * Polytope_get_a(p, i, x);}
      if(sum > Polytope_get_b(p, i)) {return false;}
   }
   return true; // passed all inequalities
}


void Polytope_intersect(const Polytope* p, const FT* x, const FT* d, FT* t0, FT* t1) {
   const int n = p->n;
   const int m = p->m;
   
   FT t00 = FT_MAX;// tmp variables for t0, t1
   FT t11 = FT_MIN;

   for(int i=0; i<m; i++) {
      // check orientation of plane vs direction of line:
      //   if orthogonal (d*ai = 0), then no intersection
      //   if >0, then same direction -> t0
      //   if <0, then opp  direction -> t1
      const FT* ai = Polytope_get_aV(p,i);
      const FT b = Polytope_get_b(p, i);
      const FT dai = dotProduct(d,ai,n);
      
      printf("dai: %f\n",dai);

      if(abs(dai) <= FT_EPS) {continue;} // orthogonal

      // find intersections y of line with all planes:
      //   y = x + d*t
      //   ai*y = b
      //   
      //   t = (b - ai*x)/(d*ai)
      
      FT t = (b - dotProduct(ai,x,n)) / dai;
      
      if(dai > 0.0) {
         t00 = (t00<t)?t00:t; // min
      } else {
         t11 = (t11<t)?t11:t; // max
      }
   }
   
   // return:
   *t0 = t00;
   *t1 = t11;
}
