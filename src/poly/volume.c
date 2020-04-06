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

Body_T Polytope_T = {
	.print = Polytope_print,
        .free = Polytope_free,
        .inside = Polytope_inside_ref,
        .intersect = Polytope_intersect_ref,
        .intersectCoord = NULL,
};
Body_T Sphere_T = {
        .print = Sphere_print,
	.free = Sphere_free,
	.inside = Sphere_inside_ref,
	.intersect = Sphere_intersect_ref,
	.intersectCoord = NULL,
};

Polytope* Polytope_new(int n, int m) {
   Polytope* o = (Polytope*) malloc(sizeof(Polytope));
   o->n = n;
   o->m = m;
   const int ftPerVec = 32/sizeof(FT); // number FT per Vector
   o->line = (n+1 + ftPerVec-1)/ftPerVec * ftPerVec;
   //o->data = (FT*) malloc(sizeof(FT)*(n+1)*m);
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



FT volumeEstimateNormalizedBody(const int n, const FT r0, const FT r1, const Polytope* body) {
   //
   // Ideas:
   //  try unit vector
   //  	 -> simplification in computation
   //  	 -> store dotProd for each ineq.
   //  try random direction
   //    -> require random from normal distr.
   //
   const int step_size = 100000; // number of points sampled
   const int walk_size = 10; // number of steps for walk
   
   // init x:
   FT* x = (FT*) malloc(sizeof(FT)*n);// sample point x
   for(int j=0;j<n;j++) {x[j]=0.0;}// origin

   FT* d = (FT*) malloc(sizeof(FT)*n); // vector for random direction

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
         for(int w=0;w<walk_size;w++) { // take some random steps for x
            int dd = prng_get_random_int_in_range(0,n-1); // pick random dimension
            for(int j=0;j<n;j++) {d[j] = ((j==dd)?1.0:0);}
            
            FT t0,t1, bt0,bt1;
            Polytope_T.intersect(body, x, d, &t0, &t1);
            Ball_intersect(n, rk, x, d, &bt0, &bt1);
            
            // ensure do not walk outside of outer ball:
            t0 = (t0>bt0)?t0:bt0; // max
            t1 = (t1<bt1)?t1:bt1; // min

            //printf("%f %f %f %f\n",bt0,bt1,t0,t1);

            FT t = prng_get_random_double_in_range(t0,t1);
            for(int j=0;j<n;j++) {x[j] += d[j]*t;}
         }
         
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
   }

   return volume;
}



FT xyz_f1(const Polytope* body, const FT r, const int n) {return body->n;}
FT xyz_f2(const Polytope* body, const FT r, const int n) {return n;}
xyz_f_t xyz_f = xyz_f1;
