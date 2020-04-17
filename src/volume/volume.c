#include "volume.h"


FT dotProduct(const FT* u, const FT* v, const int n) {
   FT sum = 0.0;
   for(int i=0; i<n; i++) {sum+= u[i]*v[i];}
   return sum;
}

void Ball_intersect(const int n, const FT r, const FT* x, const FT* d, FT* t0, FT* t1) {
   // y = x + d*t
   // y^2 = r^2
   //
   // x^2 - r^2 + 2x*d*t + d^2*t^2 = 0
   
   FT x2 = dotProduct(x,x,n);
   FT d2 = dotProduct(d,d,n);
   FT xd = dotProduct(x,d,n);

   FT a = d2;
   FT ainv = 1.0 / a;
   FT b = 2.0*xd;
   FT c = x2 - r*r;

   FT detSqrt = sqrt(b*b - 4.0*a*c);
   
   *t1 = (-b + detSqrt) * 0.5 * ainv;
   *t0 = (-b - detSqrt) * 0.5 * ainv;
}

void Ball_intersectCoord(const int n, const FT r, const FT* x, const int d, FT* t0, FT* t1) {
   FT x2 = dotProduct(x,x,n);
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

Body_T Polytope_T = {
	.print = Polytope_print,
        .free = Polytope_free,
        .inside = Polytope_inside_ref,
        .intersect = Polytope_intersect_ref,
        .intersectCoord = Polytope_intersectCoord_ref,
	.cacheAlloc = Polytope_cacheAlloc_ref,
	.cacheReset = Polytope_cacheReset_ref,
	.cacheUpdateCoord = Polytope_cacheUpdateCoord_ref,
        .shallowCutOracle = Polytope_shallowCutOracle_ref,
};
Body_T Sphere_T = {
        .print = Sphere_print,
	.free = Sphere_free,
	.inside = Sphere_inside_ref,
	.intersect = Sphere_intersect_ref,
	.intersectCoord = Sphere_intersectCoord_ref,
	.cacheAlloc = Sphere_cacheAlloc_ref,
	.cacheReset = Sphere_cacheReset_ref,
	.cacheUpdateCoord = Sphere_cacheUpdateCoord_ref,
	.shallowCutOracle = NULL,
};
Body_T Ellipsoid_T = {
        .print = Ellipsoid_print,
	.free = Ellipsoid_free,
	.inside = Ellipsoid_inside_ref,
	.intersect = Ellipsoid_intersect_ref,
	.intersectCoord = Ellipsoid_intersectCoord_ref,
	.cacheAlloc = Ellipsoid_cacheAlloc_ref,
	.cacheReset = Ellipsoid_cacheReset_ref,
	.cacheUpdateCoord = Ellipsoid_cacheUpdateCoord_ref,
	.shallowCutOracle = Ellipsoid_shallowCutOracle_ref,
};

Polytope* Polytope_new(int n, int m) {
   Polytope* o = (Polytope*) malloc(sizeof(Polytope));
   o->n = n;
   o->m = m;
   o->line = ceil_cache(n+1,sizeof(FT)); // make sure next is also 32 alligned
   o->data = (FT*)(aligned_alloc(32, o->line*m*sizeof(FT))); // align this to 32

   return o;
}


void Polytope_free(const void* o) {
   Polytope* p = (Polytope*)o;
   free(p->data);
   free(p);
}


void Polytope_print(const void* o) {
   const Polytope* p = (Polytope*)o;
   printf("Polytope: n=%d, m=%d\n",p->n,p->m);
   for(int i=0; i<p->m; i++) {
      for(int j=0; j<p->n; j++) {
         printf(" %.3f",Polytope_get_a(p,i,j));
      }
      printf(" | %.3f\n",Polytope_get_b(p,i));
   }
}

inline FT* Polytope_get_Ai(const Polytope* p, int i) {
   return &(p->data[i * (p->line)]);
}

inline void Polytope_set_a(Polytope* p, int i, int x, FT a) {
   p->data[i * (p->line) + x] = a;
}
inline void Polytope_set_b(Polytope* p, int i, FT b) {
   p->data[i * (p->line) + p->n] = b;
}

inline FT Polytope_get_a(const Polytope* p, int i, int x) {
   return p->data[i * (p->line) + x];
}
inline FT Polytope_get_b(const Polytope* p, int i) {
   return p->data[i * (p->line) + p->n];
}

bool Polytope_inside_ref(const void* o, const FT* v) {
   const Polytope* p = (Polytope*)o;
   for(int i=0; i<p->n; i++) {
      FT sum = 0;
      for(int x=0; x<p->n; x++) { sum+= v[x] * Polytope_get_a(p, i, x);}
      if(sum > Polytope_get_b(p, i)) {return false;}
   }
   return true; // passed all inequalities
}


void Polytope_intersect_ref(const void* o, const FT* x, const FT* d, FT* t0, FT* t1) {
   const Polytope* p = (Polytope*)o;
   const int n = p->n;
   const int m = p->m;
   
   FT t00 = -FT_MAX;// tmp variables for t0, t1
   FT t11 = FT_MAX;

   for(int i=0; i<m; i++) {
      // check orientation of plane vs direction of line:
      //   if orthogonal (d*ai = 0), then no intersection
      //   if <0, then same direction -> t0
      //   if >0, then opp  direction -> t1
      const FT* ai = Polytope_get_Ai(p,i);
      const FT b = Polytope_get_b(p, i);
      const FT dai = dotProduct(d,ai,n);
      // Note: if base-vector: could just pick i'th entry!
      
      //printf("dai: %f %f\n",dai,FT_EPS);

      if(dai <= FT_EPS && -dai <= FT_EPS) {continue;} // orthogonal

      // find intersections y of line with all planes:
      //   y = x + d*t
      //   ai*y = b
      //   
      //   t = (b - ai*x)/(d*ai)
      
      FT t = (b - dotProduct(ai,x,n)) / dai;
      //printf("t: %f\n",t);
      
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

void Polytope_intersectCoord_ref(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache) {
   const Polytope* p = (Polytope*)o;
   const int n = p->n;
   const int m = p->m;
   FT* Aix = (FT*)cache;
   
   FT t00 = -FT_MAX;// tmp variables for t0, t1
   FT t11 = FT_MAX;

   for(int i=0; i<m; i++) {
      const FT* ai = Polytope_get_Ai(p,i);
      const FT b = Polytope_get_b(p, i);
      const FT dai = ai[d]; // dot product with unit vector dim d
      
      if(dai <= FT_EPS && -dai <= FT_EPS) {continue;} // orthogonal
      
      const FT aix = dotProduct(ai,x,n);
      assert(aix == Aix[i] && "Cache must be accurate!");
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

void Polytope_intersectCoord_cached_ref(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache) {
   const Polytope* p = (Polytope*)o;
   const int n = p->n;
   const int m = p->m;
   FT* Aix = (FT*)cache;
   
   FT t00 = -FT_MAX;// tmp variables for t0, t1
   FT t11 = FT_MAX;

   for(int i=0; i<m; i++) {
      const FT* ai = Polytope_get_Ai(p,i);
      const FT b = Polytope_get_b(p, i);
      const FT dai = ai[d]; // dot product with unit vector dim d
      
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



int  Polytope_cacheAlloc_ref(const void* o) {
   const Polytope* p = (Polytope*)o;
   // allocate one FT per inequality: store dot product: dot(ai,x)
   return p->m * sizeof(FT);
}
void Polytope_cacheReset_ref(const void* o, const FT* x, void* cache) {
   const Polytope* p = (Polytope*)o;
   FT* c = (FT*)cache;
   const int n = p->n;
   const int m = p->m;
   for(int i=0; i<m; i++) {
      const FT* ai = Polytope_get_Ai(p,i);
      c[i] = dotProduct(ai,x,n);
   }
}

void Polytope_cacheUpdateCoord_ref(const void* o, const int d, const FT dx, void* cache) {
   const Polytope* p = (Polytope*)o;
   const int m = p->m;
   FT* c = (FT*)cache;
   for(int i=0; i<m; i++) {
      c[i] += dx * Polytope_get_a(p,i,d);
   } 
}

bool Polytope_shallowCutOracle_ref(const void* o, const Ellipsoid* e, FT* v, FT* c) {
   const Polytope* p = (Polytope*)o;
   const int n = p->n;
   const int m = p->m;
   
   FT Ax[m];

   // check if center of Ellisoid e = (X,x) is in Polytope:
   // for all i, check if:
   //    Ai * x <= bi
   for(int i=0;i<m;i++) {
      FT* Ai = Polytope_get_Ai(p,i);
      FT bi = Polytope_get_b(p,i);
      FT* x = e->a;
      Ax[i] = dotProduct(Ai,x,n);
      if(Ax[i] > bi) { // found one -> return (Ai, bi)
         for(int j=0;j<n;j++) {v[j] = Ai[j];}
	 *c = bi;
         return true;
      }
   }
   

   // check if inner Ellipsoid e = ( (2n)^-2 * T.inverse(), x) is in Polytope:
   // for all i, check if:
   //   AiT * T * Ai <= (bi - AiT * x)^2 * (2n)^2
   const FT twon2 = 4.0*n*n;
   for(int i=0;i<m;i++) {
      FT* Ai = Polytope_get_Ai(p,i);
      FT bi = Polytope_get_b(p,i);
      
      FT AitTAi = 0; // could be useful to cache...
      for(int j=0;j<n;j++) {
         FT* Tj = Ellipsoid_get_Ti(e,j);
	 FT TjAi = dotProduct(Tj,Ai,n);
         AitTAi += Ai[j] * TjAi;
      }
      
      FT diff = bi - Ax[i];
      if(AitTAi > diff*diff*twon2) { // found one -> return (Ai, bi)
         for(int j=0;j<n;j++) {v[j] = Ai[j];}
	 *c = bi;
         return true;
      }
   }
   
   // no half-plane violated inner ellipse
   return false;
}

Sphere* Sphere_new(int n, FT r, const FT* c) {
   Sphere* o = (Sphere*) malloc(sizeof(Sphere));
   o->n = n;
   o->r = r;
   o->center = (FT*)(aligned_alloc(32, n*sizeof(FT))); // align this to 32
   for(int i=0; i<n; i++) {o->center[i] = c[i];}
   return o;
}

void Sphere_free(const void* o) {
   Sphere* s = (Sphere*)o;
   free(s->center);
   free(s);
}

void Sphere_print(const void* o) {
   const Sphere* s = (Sphere*)o;
   printf("Sphere: n=%d, r=%.3f, c=[",s->n,s->r);
   for(int i=0; i<s->n; i++) {
      printf(" %.3f",s->center[i]);
   }
   printf("]\n");
}

bool Sphere_inside_ref(const void* o, const FT* v) {
   const Sphere* s = (Sphere*)o;
   FT d2 = 0.0;
   for(int i=0; i<s->n; i++) { FT d = s->center[i] - v[i]; d2 += d*d;}
   return d2 <= s->r*s->r;
}

void Sphere_intersect_ref(const void* o, const FT* x, const FT* d, FT* t0, FT* t1) {
   const Sphere* s = (Sphere*)o;
   const int n = s->n;
   FT diff[n]; // probably a terrible idea, besides not vector alligned!
   for(int i=0;i<n;i++) {diff[i] = x[i] - s->center[i];}
   Ball_intersect(n, s->r, diff, d, t0,t1);
}

void Sphere_intersectCoord_ref(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache) {
   const Sphere* s = (Sphere*)o;
   const int n = s->n;
   FT diff[n]; // probably a terrible idea, besides not vector alligned!
   for(int i=0;i<n;i++) {diff[i] = x[i] - s->center[i];}
   Ball_intersectCoord(n, s->r, diff, d, t0,t1);
}

int Sphere_cacheAlloc_ref(const void* o) {
   return 0; // no cache
}
void Sphere_cacheReset_ref(const void* o, const FT* x, void* cache) {
   // no cache
}
void Sphere_cacheUpdateCoord_ref(const void* o, const int d, const FT dx, void* cache) {
   // no cache
}

Ellipsoid* Ellipsoid_new(int n) {
   Ellipsoid* e = (Ellipsoid*) malloc(sizeof(Ellipsoid));
   e->n = n;
   e->line = ceil_cache(n,sizeof(FT)); // make sure next is also 32 alligned
   e->A = (FT*)(aligned_alloc(32, e->line*(n+1)*sizeof(FT))); // align this to 32
   e->a = e->A + e->line * n;
   e->T = NULL;
   for(int i=0; i<n; i++) {
      for(int j=0; j<n; j++) {
	 e->A[i*e->line + j] = (i==j)?1:0;
      }
      e->a[i] = 0;
   }
   return e;
}

Ellipsoid* Ellipsoid_new_with_T(int n) {
   Ellipsoid* e = (Ellipsoid*) malloc(sizeof(Ellipsoid));
   e->n = n;
   e->line = ceil_cache(n,sizeof(FT)); // make sure next is also 32 alligned
   e->A = (FT*)(aligned_alloc(32, e->line*(2*n+1)*sizeof(FT))); // align this to 32
   e->a = e->A + e->line * n;
   e->T = e->a + e->line;
   for(int i=0; i<n; i++) {
      for(int j=0; j<n; j++) {
	 e->A[i*e->line + j] = (i==j)?1:0;
	 e->T[i*e->line + j] = (i==j)?1:0;
      }
      e->a[i] = 0;
   }
   return e;
}


FT* Ellipsoid_get_Ai(const Ellipsoid* e, int i) {
   return e->A + i*e->line;
}

FT* Ellipsoid_get_Ti(const Ellipsoid* e, int i) {
   assert(e->T && "T must be allocated");
   return e->T + i*e->line;
}

void Ellipsoid_free(const void* o) {
   Ellipsoid* e = (Ellipsoid*)o;
   free(e->A);// includes a (and T)
   free(e);
}

void Ellipsoid_print(const void* o) {
   Ellipsoid* e = (Ellipsoid*)o;
   const int n = e->n;
   printf("Ellipsoid: n=%d\n",e->n);
   printf("A:\n");
   for(int i=0;i<n;i++) {
      const FT* Ai = Ellipsoid_get_Ai(e,i);
      for(int j=0;j<n;j++) {
         printf("%.12f ",Ai[j]);
      }
      printf("\n");
   } 
   if(e->T){
      printf("T:\n");
      for(int i=0;i<n;i++) {
         const FT* Ti = Ellipsoid_get_Ti(e,i);
         for(int j=0;j<n;j++) {
            printf("%.12f ",Ti[j]);
         }
         printf("\n");
      }
   }
   printf("a:\n");
   for(int j=0;j<n;j++) {
      printf("%.12f ",e->a[j]);
   }
   printf("\n");
 
}

bool Ellipsoid_inside_ref(const void* o, const FT* v) {
   Ellipsoid* e = (Ellipsoid*)o;
   return Ellipsoid_eval(e, v) <= 1.0;
}

void Ellipsoid_intersect_ref(const void* o, const FT* x, const FT* d, FT* t0, FT* t1) {
   // (y-a)T * A * (y-a) = 1
   // y = x + t*d
   //
   // z = x-a
   // t^2 * (dT * A * d) + t * 2(dT * A * z) + (zT * A * z) - 1 = 0
   //
   // Simple quadratic root problem
   // t^2 * a + t * b + c = 0
   
   Ellipsoid* e = (Ellipsoid*)o;
   const int n = e->n;

   FT a = 0;
   FT b = 0;
   FT c = -1.0;
   
   // do multiplications same as in eval.
   for(int i=0;i<n;i++) {
      const FT* Ai = Ellipsoid_get_Ai(e,i);
      FT Ad = 0;
      FT Az = 0;
      for(int j=0; j<n; j++) {
         Ad += Ai[j] * d[j];
         Az += Ai[j] * (x[j] - e->a[j]);
      }
      a += d[i] * Ad;
      b += d[i] * Az;
      c += (x[i] - e->a[i]) * Az;
   }
   b *= 2.0;

   // find t:
   const FT det = b*b - 4.0*a*c;
   assert(det >= 0);
   const FT sqrtDet = sqrt(det);
   const FT aInv = 0.5/a;

   *t0 = (-b - sqrtDet) * aInv;
   *t1 = (-b + sqrtDet) * aInv;
}

void Ellipsoid_intersectCoord_ref(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache) {
   Ellipsoid* e = (Ellipsoid*)o;
   assert(false && "not implemented!");
}

int  Ellipsoid_cacheAlloc_ref(const void* o) {
   Ellipsoid* e = (Ellipsoid*)o;
   assert(false && "not implemented!");
}

void Ellipsoid_cacheReset_ref(const void* o, const FT* x, void* cache) {
   Ellipsoid* e = (Ellipsoid*)o;
   assert(false && "not implemented!");
}

void Ellipsoid_cacheUpdateCoord_ref(const void* o, const int d, const FT dx, void* cache) {
   Ellipsoid* e = (Ellipsoid*)o;
   assert(false && "not implemented!");
}

bool Ellipsoid_shallowCutOracle_ref(const void* o, const Ellipsoid* e, FT* v, FT* c) {
   Ellipsoid* this = (Ellipsoid*)o;
   const int n = this->n;

   // case 1: center of cage (e->a) is outside of body o:
   if(!Ellipsoid_T.inside(this, e->a)) {
      //printf("not inside ellipsoid!\n");
      
      Ellipsoid_normal(this, e->a, v);
      *c = dotProduct(v,e->a, n);

      return true;
   }

   // run minimization to obtain a point where to cut:
   FT* x0 = (FT*)(aligned_alloc(32, n*sizeof(FT))); // align this to 32
   FT* x1 = (FT*)(aligned_alloc(32, n*sizeof(FT))); // align this to 32

   for(int i=0;i<n;i++) {x0[i]=prng_get_random_double_normal();}// random direction
   for(int i=0;i<n;i++) {x1[i]=this->a[i] + x0[i];}
   for(int i=0;i<n;i++) {x0[i]=this->a[i] + x0[i];}
   FT beta2 = 1.0 / (4*n*n);
   Ellipsoid_minimize(this,1, e, x0);
   Ellipsoid_minimize(this,1, e, x1);
   
   FT eval0 = Ellipsoid_eval(e,x0); 
   FT eval1 = Ellipsoid_eval(e,x1);
   
   //printf("eval: %f %f vs %f\n",eval0,eval1,beta2);
   if(eval0 > beta2 && eval1 > beta2) { return false; } // both local minima too far out
   
   // choose better x:
   FT* x = x0;
   if(eval1 < eval0) {x = x1;}
   
   Ellipsoid_normal(this, x, v);
   *c = dotProduct(v,x, n);
}

FT Ellipsoid_eval(const Ellipsoid* e, const FT* x) {
   // (x-a)T * A * (x-a)
   const int n = e->n;
   FT sum = 0;
   for(int i=0;i<n;i++) {
      const FT* Ai = Ellipsoid_get_Ai(e,i);
      FT Axa = 0;
      for(int j=0; j<n; j++) {
         Axa += Ai[j] * (x[j] - e->a[j]);
      }
      sum += (x[i] - e->a[i]) * Axa;
   }
   return sum;
}

void Ellipsoid_normal(const Ellipsoid* e, const FT* x, FT* normal) {
   int n = e->n;
   for(int i=0;i<n;i++) {
      const FT* Ai = Ellipsoid_get_Ai(e,i);
      FT Axa = 0;
      for(int j=0; j<n; j++) {
         Axa += Ai[j] * (x[j] - e->a[j]);
      }
      normal[i] = 2.0*Axa;
   }
}

void Ellipsoid_project(const Ellipsoid* e, const FT eFac, FT* x) {
   // internal.
   // push x back on surface of e
   // 
   // for now just pull to center. Could try with normal also...?
   int n = e->n;
   FT eval = Ellipsoid_eval(e,x);
   FT scale = sqrt(eFac / eval);
   for(int i=0; i<n; i++) { x[i] = e->a[i] + (x[i]-e->a[i]) * scale;}
}

void Ellipsoid_minimize(const Ellipsoid* e, const FT eFac, const Ellipsoid* f, FT* x){
   // Argument why we cannot have more than 2 strict local minima
   // 
   // First, reduce to 2d problem
   // Assume you have 3 strict local minima in n-dim problem,
   // take intersection with ellipsoid and 2d plane of the three points, get 2d ellipse.
   // Now you have 3 strict local minima on a 2d ellipse for a 2d ellipse cost function.
   // 
   // Why you cannot have 3 strict local minima in 2d case:
   // 
   // scale cost function to unit circle.
   // for solution points, the normals must be parallel
   // B * (x - b) = lambda * I * x
   // Prove this only holds for 4 points: 2 local minima and 2 local maxima
   // or for all points, all min=max
   //
   // B * (x - b) = lambda * I * x
   //
   // (B - lambda * I) * x = B*b
   //
   // Almost looks like eigenvector thing... but not with =0
   // We are not sure how to make this formal...

   const int n = e->n;
   
   // can we alloc before somehow?
   FT* nE = (FT*)(aligned_alloc(32, n*sizeof(FT))); // align this to 32
   FT* nF = (FT*)(aligned_alloc(32, n*sizeof(FT))); // align this to 32
   FT* nP = (FT*)(aligned_alloc(32, n*sizeof(FT))); // align this to 32
   
   // for debugging:
   FT* tmp = (FT*)(aligned_alloc(32, n*sizeof(FT))); // align this to 32

   Ellipsoid_project(e,eFac,x);
   
   int count = 0;
   FT dot = FT_MAX; // step size
   FT eval = FT_MAX; // step size
   FT beta = 1.0;
   int lastBeta = 0;
   do{
      Ellipsoid_normal(e, x, nE);
      FT nE2 = dotProduct(nE,nE, n);
      Ellipsoid_normal(f, x, nF);
      
      FT proj = dotProduct(nE,nF,n)/nE2;
      // project gradient in f on plane (normal to e)
      for(int i=0; i<n; i++) {
         nP[i] = nF[i] - nE[i]*proj;
      }
      FT lastDot = dot;
      dot = dotProduct(nP,nP, n);
      
      // debug output
      FT lastEval = eval;
      eval = Ellipsoid_eval(f,x);
      //printf("dot %.12f eval %.12f\n",dot,eval);
      // adjust step rate beta
      if(lastEval < eval) {
         beta *= 0.5;
	 //printf("beta down! %.20f\n",beta);
	 lastBeta = 0;
      }
      if(lastBeta++ > 10) {
	 beta*=1.2;
	 //printf("beta up! %.20f\n",beta);
	 lastBeta = 0;
      }
      
      for(int i=0; i<n; i++) {tmp[i] = x[i];}// debug
      for(int i=0; i<n; i++) {
         x[i] -= beta * nP[i]; // step = c - n (n * c)
      }

      Ellipsoid_project(e,eFac,x);// pull back on ellipsoid

      //for(int i=0; i<n; i++) {tmp[i] -= x[i];}// debug
      //FT dotTmp = dotProduct(tmp,tmp,n);
      //printf("moved: %.20f\n",dotTmp);
      //FT evalX = Ellipsoid_eval(e,x);
      //printf("evalX: %.20f\n",evalX);
      
      if(++count >= 100*n) {
         printf("Warning: taking too many steps (%d), abort minimization now!\n",count);
	 break;
      }
   } while(dot > 0.00000001);
   
   printf("steps taken: %d\n",count);
   
   //for(int i=0; i<n; i++) {printf(" %.12f", x[i]);}// debug
   //printf("\n");

   free(nE);
   free(nF);
   free(nP);
}


void preprocess_ref(const int n, const int bcount, const void** body_in, void** body_out, const Body_T** type, FT *det) {
   // 1. init_ellipsoid:
   //     idea: origin at 0, radius determined by min of all bodies
   printf("init_ellipsoid\n");
   
   FT R2 = 1e6; // TODO - ask sub bodies for radius, take min
   
   Ellipsoid* e = Ellipsoid_new_with_T(n); // origin zero
   for(int i=0; i<n; i++) {
      FT* Ai = Ellipsoid_get_Ai(e,i);
      FT* Ti = Ellipsoid_get_Ti(e,i);
      Ai[i] = 1.0/R2; // sphere with 
      Ti[i] = R2; // sphere with 
   }
   Ellipsoid_T.print(e);
   
   // 2. Cut steps
   printf("cut steps\n");
   
   const FT beta_r = 2*n;
   const FT beta = 1.0/beta_r;

   const FT ro = (1.0 - n*beta) / (n+1.0);// c2
   const FT tow = 2.0 * ro / (1.0 - beta); // c4
   const FT onemnb = (1.0 - n*beta);
   const FT zs = (2.0*n*n + onemnb*onemnb) * (1.0 - beta*beta) / (2.0*n*n-2.0); // zeta*sigma 

   FT c; // plane for oracle: normal*x <= x
   FT* v = (FT*)aligned_alloc(32, n*sizeof(FT)); // align this to 32

   int step = 0;
   while(true) {
      step++;
      //printf("cut step %d\n",step);
      
      bool doCut = false;
      for(int b=0; b<bcount; b++) {
         doCut = type[b]->shallowCutOracle(body_in[b], e, v, &c);
	 if(doCut) {break;}
      }

      if(!doCut) {break;} // we are done!

      //printf("cut!\n");
      
      // calculate: T * v and vt * T * v
      FT Tv[n];
      FT vtTv = 0;
      for(int i=0;i<n;i++) {
         const FT* Ti = Ellipsoid_get_Ti(e,i);
	 Tv[i] = dotProduct(Ti,v,n);
	 vtTv += v[i] * Tv[i];
      }
      
      // update a:
      FT* a = e->a;
      FT fac = ro/sqrt(vtTv);
      for(int i=0;i<n;i++) {
         a[i] -= fac * Tv[i];
      }

      // update T:
      // T' = zs * (T - fac2*Tv*Tvt)
      FT fac2 = tow / vtTv;
      for(int i=0;i<n;i++) {
         FT* Ti = Ellipsoid_get_Ti(e,i);
         for(int j=0;j<n;j++) {
	    Ti[j] = zs * (Ti[j] - fac2*Tv[i]*Tv[j]);
	 }
      }

      // update A:
      // above, we made a rank-1 update for T.
      // We can update A (=T.inverse()) using the Sherman-Morisson formula:
      // (T + u*vt).inverse() = T.inverse() - T.inverse() * u * vt * T.inverse() / (1 + vt * T.inverse() * u)
      // A' = (T + u*vt).inverse() = A - A * u * vt * A / (1 + vt * A * u)
      //
      // So we will compute:
      // A' = A - fac2 * A  * Tv * Tvt * A / (1 + fac2 * Tvt * A * Tv)
      
      FT ATv[n];
      FT TvtATv = 0;
      for(int i=0;i<n;i++) {
         FT* Ai = Ellipsoid_get_Ai(e,i);
         ATv[i] = 0;
         for(int j=0;j<n;j++) {// dotProduct
	    ATv[i] += Ai[j] * Tv[j];
	 }

	 TvtATv += Tv[i] * ATv[i];
      }
      FT div = 1.0 / (1 + fac2*TvtATv);

      for(int i=0;i<n;i++) {
         FT* Ai = Ellipsoid_get_Ai(e,i);
         for(int j=0;j<n;j++) {
	    Ai[j] = (Ai[j] - fac2 * ATv[i]*ATv[j]) * div;
	 }
      } 
   }
   Ellipsoid_T.print(e);
   printf("took %d steps.\n",step);

   // 3. Transformation
   
   // outer ellipsoid (cage):
   // (x-a)T * T.inverse() * (x-a) <= 1
   // 
   // We want to transform our space, such that this ellipsoid becomes a ball:
   // (y-a') * I * (y - a') <=1
   // and a' == 0
   // 
   // Cholesky:
   // T = L * LT
   // 
   // How do we get there:
   //
   // (x-a)T * (L * LT).inverse() * (x-a) <= 1
   // (x-a)T * LT.inverse() * L.inverse() * (x-a) <= 1
   // (L.inverse()*x-L.inverse()*a)T * (L.inverse()*x-L.inverse()*a) <= 1
   // 
   // y = L.invere()*x - L.inverse()*a
   // a' = 0
   // x = L * y + a
   //
   // Polytope:
   // A * x <= b
   // A * (L * y + a) <= b
   // A * L * y <= b - A * a
   // A' = A * L
   // b' = b - A * a
   // A * L * y <= b'
   //
   // Paper then grows everything by 1/beta, such that unit ball = inner ellipsoid
   // hence b'' = b' / beta
   //       A'' = A'
   // 
   // Ellipsoid:
   // (x-b)T * B * (x-b) <= 1
   // (L * y + a - b)T * B * (L * y + a - b) <= 1
   // 
   // goal:
   // (y - b')T * B' * (y - b') <= 1
   //
   // B' = LT * B * L
   // b' = L.inverse() * (a-b)
   //
   // grow by 1/beta:
   // (y - b')T * B' * (y - b') <= 1/beta
   // B'' = B' * beta 
   // b'' = b'

   printf("Transform\n");
   //FT *L = (FT *) calloc(n*n, sizeof(FT));
   //int err = cholesky(T, L, n);
   //
   //if (err > 0){
   //   printf("The input polytope is degenerate or non-existant and the volume is 0.\n");
   //   exit(1);		
   //}
   
   
   //// b = beta_r * (b - A * ori);
   //for (int i = 0; i < m; i++){
   //   Polytope_set_b(*Q, i, beta_r * distance[i]);
   //}
   //// A = A * Trans;
   //for (int i = 0; i < m; i++){
   //   for (int j = 0; j < n; j++){
   //      FT sum = 0;
   //      for (int k = 0; k < n; k++){
   //         sum += Polytope_get_a(P, i, k) * Trans[k*n + j];
   //      }
   //      Polytope_set_a(*Q, i, j, sum);
   //   }
   //}


   //*det = 1;
   //for (int i = 0; i < n; i++){
   //   *det *= Trans[i*n+i];
   //}
   //*det /= pow(beta_r, n);
}


int step_size = 100000;
int walk_size = 10;

walk_f_t walk_f = walk_ref;

void walk_ref(const int n, const FT rk, int bcount, const void** body, const Body_T** type, FT* x, FT* d, void** cache) {
   const int ws = walk_size; // number of steps for walk
   
   for(int w=0;w<ws;w++) { // take some random steps for x
      // set d to random direction vector, not normalized
      for(int j=0;j<n;j++) {d[j] = prng_get_random_double_normal();}
      
      FT t0,t1;
      // ensure do not walk outside of outer ball:
      Ball_intersect(n, rk, x, d, &t0, &t1);
      
      for(int c=0;c<bcount;c++) {
         FT bt0, bt1;
         type[c]->intersect(body[c], x, d, &bt0, &bt1);
         t0 = (t0>bt0)?t0:bt0; // max
         t1 = (t1<bt1)?t1:bt1; // min
      }
   
      FT t = prng_get_random_double_in_range(t0,t1);
      for(int j=0;j<n;j++) {x[j] += d[j]*t;}
   }
}

void walkCoord_ref(const int n, const FT rk, int bcount, const void** body, const Body_T** type, FT* x, FT* d, void** cache) {
   const int ws = walk_size; // number of steps for walk
   
   for(int w=0;w<ws;w++) { // take some random steps for x
      int dd = prng_get_random_int_in_range(0,n-1); // pick random dimension
      
      FT t0,t1;
      // ensure do not walk outside of outer ball:
      Ball_intersectCoord(n, rk, x, dd, &t0, &t1);
      
      for(int c=0;c<bcount;c++) {
         FT bt0, bt1;
         type[c]->intersectCoord(body[c], x, dd, &bt0, &bt1, cache[c]);
         t0 = (t0>bt0)?t0:bt0; // max
         t1 = (t1<bt1)?t1:bt1; // min
      }
   
      FT t = prng_get_random_double_in_range(t0,t1);
      x[dd] += t;
      for(int c=0;c<bcount;c++) {
         type[c]->cacheUpdateCoord(body[c], dd, t, cache[c]);
      }

   }
}



FT volume_ref(const int n, const FT r0, const FT r1, int bcount, const void** body, const Body_T** type) {
   //
   // Ideas:
   //  try unit vector
   //  	 -> simplification in computation
   //  	 -> store dotProd for each ineq.
   //  try random direction
   //    -> require random from normal distr.
   //
   //  cache values for bodies:
   //    -> body has function saying how much memory.
   //    -> body has function to initialize it for an x
   //    -> how update? via function or direct?
   //    -> allocate in beginning, give offsets for each body
   //
   
   // init x:
   FT* x = (FT*) malloc(sizeof(FT)*n);// sample point x
   for(int j=0;j<n;j++) {x[j]=0.0;}// origin

   FT* d = (FT*) malloc(sizeof(FT)*n); // vector for random direction

   // set up cache:
   int cache_size = 0;
   void* cache[bcount];
   for(int c=0;c<bcount;c++) {
      cache_size += type[c]->cacheAlloc(body[c]);
   }
   void* cache_base = aligned_alloc(32, cache_size); // align this to 32
   cache_size = 0;
   for(int c=0;c<bcount;c++) {
      cache[c] = cache_base + cache_size;
      cache_size += type[c]->cacheAlloc(body[c]);

      type[c]->cacheReset(body[c],x,cache[c]);
   }

   
   const int l = ceil(n*log(r1/r0) / log(2.0));
   printf("steps: %d\n",l);
   int t[l+1];// counts how many were thrown into Bi
   for(int i=0;i<=l;i++){t[i]=0;}
   
   // volume up to current step
   // start with B(0,r0)
   // multiply with estimated factor each round
   FT volume = Ball_volume(n, r0);
   
   // Idea:
   //   last round must end in rk = r0/stepFac
   //   first round must start with rk >= r1
   const FT stepFac = pow(2,-1.0/(FT)n);
   
   FT rk = r0*pow(stepFac,-l);
   int count = 0;
   for(int k=l;k>0;k--,rk*=stepFac) { // for each Bk
      FT kk = log(rk/r0)/(-log(stepFac));
      printf("k: %d rk: %f kk: %f step: %f\n",k,rk,kk,log(stepFac));

      //{
      //   FT rtest = (rk*rk)*0.99;
      //   FT m = log(rtest/(r0*r0))*0.5/(-log(stepFac));
      //   printf("m: %f\n",m);
      //}
      //{
      //   FT rtest = (r0*r0)*0.99;
      //   FT m = log(rtest/(r0*r0))*0.5/(-log(stepFac));
      //   printf("m: %f\n",m);
      //}
      for(int i=count; i<step_size; i++) { // sample required amount of points
         walk_f(n, rk, bcount, body, type, x, d, (void**)(&cache));
        
         // find right Bm:
         const FT x2 = dotProduct(x,x,n); // normalized radius
         const FT mmm = log(x2/(r0*r0))*0.5/(-log(stepFac));
         const int mm = ceil(mmm);
	 const int m = (mm>0)?mm:0; // find index of balls

	 //printf("k %d  m %d\n",k,m);
         assert(m <= k);
         t[m]++;
      }

      // update count:
      count = 0;
      for(int i=0;i<k;i++){count+=t[i];}
      // all that fell into lower balls

      FT ak = (FT)step_size / (FT)count;
      volume *= ak;

      printf("count: %d, volume: %f\n",count,volume);

      // x = stepFac * x   -> guarantee that in next smaller ball
      for(int j=0;j<n;j++) {x[j] *= stepFac;}
      // may have to set to zero if initializing simpler for that
      for(int c=0;c<bcount;c++) {
         type[c]->cacheReset(body[c],x,cache[c]);
      }
   }

   free(x);
   free(d);
   free(cache_base);

   return volume;
}



FT xyz_f1(const Polytope* body, const FT r, const int n) {return body->n;}
FT xyz_f2(const Polytope* body, const FT r, const int n) {return n;}
xyz_f_t xyz_f = xyz_f1;




Matrix* Matrix_new(int n, int m) {
   Matrix* o = (Matrix*) malloc(sizeof(Matrix));
   o->n = n;
   o->m = m;
   o->line = ceil_cache(n,sizeof(FT)); // make sure next is also 32 alligned
   o->data = (FT*)(aligned_alloc(32, o->line*m*sizeof(FT))); // align this to 32

   return o;
}


void Matrix_free(const void* o) {
   Matrix* p = (Matrix*)o;
   free(p->data);
   free(p);
}

inline FT* Matrix_get_row(const Matrix* m, int i) {
   return &(m->data[i * (m->line)]);
}

inline void Matrix_set(Matrix* m, int i, int x, FT a) {
   m->data[i * (m->line) + x] = a;
}
inline FT Matrix_get(const Matrix* m, int i, int x) {
   return m->data[i * (m->line) + x];
}

void Matrix_print(const void* o) {
   const Matrix* p = (Matrix*)o;
   printf("Matrix: n=%d, m=%d\n",p->n,p->m);
   for(int i=0; i<p->m; i++) {
      for(int j=0; j<p->n; j++) {
         printf(" %.3f", Matrix_get(p,i,j));
      }
      printf("\n");
   }
}
