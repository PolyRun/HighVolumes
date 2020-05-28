#include "polytopeT.h"

Body_T PolytopeT_T = {
	.print = PolytopeT_print,
        .free = PolytopeT_free,
        .clone = PolytopeT_clone,
        .inside = PolytopeT_inside_ref,
        .intersect = PolytopeT_intersect_ref,
        .intersectCoord = PolytopeT_intersectCoord_ref,
	.cacheAlloc = PolytopeT_cacheAlloc_ref,
	.cacheReset = PolytopeT_cacheReset_ref,
	.cacheUpdateCoord = PolytopeT_cacheUpdateCoord_ref,
        .shallowCutOracle = PolytopeT_shallowCutOracle_ref,
	.transform = PolytopeT_transform_ref,
        .boundingSphere = PolytopeT_bounding_ref,
        .intersectCoord4 = PolytopeT_intersectCoord4_ref,
        .intersectCoord8 = PolytopeT_intersectCoord8_ref,
	.cacheReset4 = PolytopeT_cacheReset4_ref,
	.cacheReset8 = PolytopeT_cacheReset8_ref,
	.cacheUpdateCoord4 = PolytopeT_cacheUpdateCoord4_ref,
	.cacheUpdateCoord8 = PolytopeT_cacheUpdateCoord8_ref,
};

PolytopeT* PolytopeT_new(int n, int m) {
   PolytopeT* o = (PolytopeT*) malloc(sizeof(PolytopeT));
   o->n = n;
   o->m = m;
   o->line = ceil_cache(m,sizeof(FT)); // make sure next is also 32 alligned
   int size_A = o->line*n;
   o->A = (FT*)(aligned_alloc(32, (2*size_A+o->line)*sizeof(FT))); // align this to 32
   // MB: test setting invalid elements to 0 for vectorization
   memset(o->A, 0, size_A*sizeof(FT));
   for(int i=0;i<size_A+o->line;i++) {o->A[i]=0;}
   o->b = o->A + size_A;
   o->Ainv = o->b + o->line;
   memset(o->Ainv, 0, size_A*sizeof(FT));
   for(int i=0;i<size_A;i++) {o->Ainv[i]=0;}
   return o;
}

void PolytopeT_fix_inv(const void* o) {
   PolytopeT* p = (PolytopeT*)o;
   const int n = p->n;
   const int m = p->m;
   for(int i=0; i<m; i++) {
      for(int j=0; j<n; j++) {
         FT aij = PolytopeT_get_a(p,i,j);
         if(aij == 0) {
	    PolytopeT_set_aInv(p, i,j, 0.0/0.0); // nan
	 } else {
	    PolytopeT_set_aInv(p, i,j, 1.0/aij);
	 }
      }
   }
}

void PolytopeT_free(const void* o) {
   PolytopeT* p = (PolytopeT*)o;
   free(p->A);
   free(p);
}

void* PolytopeT_clone(const void* o) {
   PolytopeT* old = (PolytopeT*)o;
   const int n = old->n;
   const int m = old->m;
   PolytopeT* p = PolytopeT_new(n,m);
   
   for(int i=0; i<m; i++) {
      for(int j=0; j<n; j++) {
         PolytopeT_set_a(p, i,j, PolytopeT_get_a(old,i,j));
         PolytopeT_set_aInv(p, i,j, PolytopeT_get_aInv(old,i,j));
      }
      PolytopeT_set_b(p, i, PolytopeT_get_b(old,i));
   }
   return p;
}

void PolytopeT_print(const void* o) {
   const PolytopeT* p = (PolytopeT*)o;
   printf("PolytopeT: n=%d, m=%d\n",p->n,p->m);
   for(int i=0; i<p->m; i++) {
      for(int j=0; j<p->n; j++) {
         printf(" %.3f",PolytopeT_get_a(p,i,j));
      }
      printf(" | %.3f\n",PolytopeT_get_b(p,i));
   }
   printf("  - inv:\n");
   for(int i=0; i<p->m; i++) {
      for(int j=0; j<p->n; j++) {
         printf(" %.3f",PolytopeT_get_aInv(p,i,j));
      }
      printf(" | %.3f\n",PolytopeT_get_b(p,i));
   }
}

bool PolytopeT_inside_ref(const void* o, const FT* v) {
   const PolytopeT* p = (PolytopeT*)o;
   for(int i=0; i<p->m; i++) {
      FT sum = 0;
      for(int x=0; x<p->n; x++) { sum+= v[x] * PolytopeT_get_a(p, i, x);}
      if(sum > PolytopeT_get_b(p, i)) {return false;}
   }
   return true; // passed all inequalities
}


void PolytopeT_intersect_ref(const void* o, const FT* x, const FT* d, FT* t0, FT* t1) {
   const PolytopeT* p = (PolytopeT*)o;
   const int n = p->n;
   const int m = p->m;
   
   FT t00 = -FT_MAX;// tmp variables for t0, t1
   FT t11 = FT_MAX;

   for(int i=0; i<m; i++) {
      // check orientation of plane vs direction of line:
      //   if orthogonal (d*ai = 0), then no intersection
      //   if <0, then same direction -> t0
      //   if >0, then opp  direction -> t1
      const FT b = PolytopeT_get_b(p, i);
      FT dai = 0;// dotProduct(d,ai,n);
      FT aix = 0;// dotProduct(x,ai,n);
      for(int j=0;j<n;j++) {
         dai+= d[j] * PolytopeT_get_a(p,i,j);
         aix+= x[j] * PolytopeT_get_a(p,i,j);
      }
      

      if(dai <= FT_EPS && -dai <= FT_EPS) {continue;} // orthogonal

      // find intersections y of line with all planes:
      //   y = x + d*t
      //   ai*y = b
      //   
      //   t = (b - ai*x)/(d*ai)
      
      FT t = (b - aix) / dai;
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

void PolytopeT_intersectCoord_ref(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache) {
   const PolytopeT* p = (PolytopeT*)o;
   const int n = p->n;
   const int m = p->m;
   FT* Aix = (FT*)cache;
   
   FT t00 = -FT_MAX;// tmp variables for t0, t1
   FT t11 = FT_MAX;
   
   for(int i=0; i<m; i++) {
      const FT b = PolytopeT_get_b(p, i);
      const FT dai = PolytopeT_get_a(p,i,d); // dot product with unit vector dim d
      
      if(dai <= FT_EPS && -dai <= FT_EPS) {continue;} // orthogonal
      
      FT aix = 0;//dotProduct(ai,x,n);
      for(int j=0;j<n;j++) {
         aix += x[j] * PolytopeT_get_a(p,i,j);
      }
      assert(abs(aix - Aix[i]) - 0.0000001 && "Cache must be accurate!");
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

int  PolytopeT_cacheAlloc_ref(const void* o) {
   const PolytopeT* p = (PolytopeT*)o;
   // allocate one FT per inequality: store dot product: dot(ai,x)
   return p->m * sizeof(FT);
}
void PolytopeT_cacheReset_ref(const void* o, const FT* x, void* cache) {
   const PolytopeT* p = (PolytopeT*)o;
   FT* c = (FT*)cache;
   const int n = p->n;
   const int m = p->m;
   for(int i=0; i<m; i++) {
      FT dot = 0;
      for(int j=0;j<n;j++) {
         dot += x[j] * PolytopeT_get_a(p,i,j);
      }
      c[i] = dot;
   }
}

void PolytopeT_cacheUpdateCoord_ref(const void* o, const int d, const FT dx, void* cache) {
   const PolytopeT* p = (PolytopeT*)o;
   const int m = p->m;
   FT* c = (FT*)cache;
   for(int i=0; i<m; i++) {
      c[i] += dx * PolytopeT_get_a(p,i,d);
   } 
}

bool PolytopeT_shallowCutOracle_ref(const void* o, const Ellipsoid* e, FT* v, FT* c) {
   const PolytopeT* p = (PolytopeT*)o;
   const int n = p->n;
   const int m = p->m;
   
   FT Ax[m];

   // check if center of Ellisoid e = (X,x) is in PolytopeT:
   // for all i, check if:
   //    Ai * x <= bi
   
   int i0 = prng_get_random_int_in_range(0,m-1);// just an experiment to see if it helps balance things
   for(int ii=0;ii<m;ii++) {
      int i = (ii+i0) % m;
      FT bi = PolytopeT_get_b(p,i);
      FT* x = e->a;
      FT Axi = 0;// dotProduct(Ai,x,n);
      for(int j=0;j<n;j++) {Axi += x[j] * PolytopeT_get_a(p,i,j);}
      Ax[i] = Axi;
      if(Ax[i] > bi) { // found one -> return (Ai, bi)
         for(int j=0;j<n;j++) {v[j] = PolytopeT_get_a(p,i,j);}
         *c = bi;
         return true;
      }
   }
   

   // check if inner Ellipsoid e = ( (2n)^-2 * T.inverse(), x) is in PolytopeT:
   // for all i, check if:
   //   AiT * T * Ai <= (bi - AiT * x)^2 * (2n)^2
   const FT twon2 = 4.0*n*n;
   int i1 = prng_get_random_int_in_range(0,m-1);//ballance experiment
   for(int ii=0;ii<m;ii++) {
      int i = (ii+i1) % m;
      FT bi = PolytopeT_get_b(p,i);
      
      FT AitTAi = 0; // could be useful to cache...
      for(int j=0;j<n;j++) {
         FT* Tj = Ellipsoid_get_Ti(e,j);
         FT TjAi = 0;// dotProduct(Tj,Ai,n);
         for(int k=0;k<n;k++) { TjAi += Tj[k] * PolytopeT_get_a(p,i,k);}
	 AitTAi += TjAi * PolytopeT_get_a(p,i,j);
      }
      
      FT diff = bi - Ax[i];
      if(AitTAi > diff*diff*twon2) { // found one -> return (Ai, bi)
         for(int j=0;j<n;j++) {v[j] = PolytopeT_get_a(p,i,j);}
         *c = bi;
         return true;
      }
   }
   
   // no half-plane violated inner ellipse
   return false;
}

void PolytopeT_transform_ref(const void* o_in, void* o_out, const Matrix* L, const FT* a, const FT beta) {
   const PolytopeT* p_in = (PolytopeT*)o_in;
   PolytopeT* p_out = (PolytopeT*)o_out;
   const int n = p_in->n;
   const int m = p_in->m;
   // computation according to explanation in preprocess_ref
   
   // b' = b - A * a
   // b'' = b' / beta
   FT beta_r = 1.0 / beta; 
   for (int i = 0; i < m; i++){
      //FT* Ai = PolytopeT_get_Ai(p_in,i);
      FT bi = PolytopeT_get_b(p_in, i);
      FT distance = bi;// - dotProduct(Ai, a, n);
      for(int j=0;j<n;j++) {distance -= a[j] * PolytopeT_get_a(p_in,i,j);}
      PolytopeT_set_b(p_out, i, beta_r * distance);
   }
   
   // A'' = A' = A * L
   for (int i = 0; i < m; i++){
      for (int j = 0; j < n; j++){
         FT sum = 0;
         for (int k = 0; k < n; k++){
            sum += PolytopeT_get_a(p_in, i, k) * Matrix_get(L, k, j);
         }
         PolytopeT_set_a(p_out, i, j, sum);
      }
   }
}
