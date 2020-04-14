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

inline FT* Polytope_get_aV(const Polytope* p, int i) {
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
      const FT* ai = Polytope_get_aV(p,i);
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
      const FT* ai = Polytope_get_aV(p,i);
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
      const FT* ai = Polytope_get_aV(p,i);
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
      const FT* ai = Polytope_get_aV(p,i);
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
   for(int i=0; i<n; i++) {
      for(int j=0; j<n; j++) {
	 e->A[i*e->line + j] = (i==j)?1:0;
      }
      e->a[i] = 0;
   }
   return e;
}

FT* Ellipsoid_get_Ai(const Ellipsoid* e, int i) {
   return e->A + i*e->line;
}

void Ellipsoid_free(const void* o) {
   Ellipsoid* e = (Ellipsoid*)o;
   free(e->A);// includes a
   free(e);
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

void Ellipsoid_project(const Ellipsoid* e, FT* x) {
   // internal.
   // push x back on surface of e
   // 
   // for now just pull to center. Could try with normal also...?
   int n = e->n;
   FT eval = Ellipsoid_eval(e,x);
   FT scale = 1.0/sqrt(eval);
   for(int i=0; i<n; i++) { x[i] = e->a[i] + (x[i]-e->a[i]) * scale;}
}

void Ellipsoid_minimize(const Ellipsoid* e, const Ellipsoid* f, FT* x){
   const int n = e->n;
   
   // can we alloc before somehow?
   FT* nE = (FT*)(aligned_alloc(32, n*sizeof(FT))); // align this to 32
   FT* nF = (FT*)(aligned_alloc(32, n*sizeof(FT))); // align this to 32
   FT* nP = (FT*)(aligned_alloc(32, n*sizeof(FT))); // align this to 32

   Ellipsoid_project(e,x);
   
   int count = 0;
   FT dot = FT_MAX; // step size
   FT beta = 0.1;
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
      
      // adjust step rate beta
      if(lastDot < dot) {beta *= 0.5; printf("beta down! %f\n",beta); lastBeta = 0;}
      if(lastBeta++ > 10) {beta*=1.2; printf("beta up! %f\n",beta); lastBeta = 0;}
      
      // debug output
      FT eval = Ellipsoid_eval(f,x);
      printf("hello %f %f\n",dot,eval);
      
      for(int i=0; i<n; i++) {
         x[i] -= beta * nP[i]; // step = c - n (n * c)
      }

      Ellipsoid_project(e,x);
   } while(count++ < 100*n && dot > 0.0000001);
   printf("steps taken: %d\n",count);
   free(nE);
   free(nF);
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
