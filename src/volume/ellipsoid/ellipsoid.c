#include "ellipsoid.h"


Body_T Ellipsoid_T = {
        .print = Ellipsoid_print,
	.free = Ellipsoid_free,
	.clone = Ellipsoid_clone_with_T,
	.inside = Ellipsoid_inside_ref,
	.intersect = Ellipsoid_intersect_ref,
	.intersectCoord = Ellipsoid_intersectCoord_ref,
	.cacheAlloc = Ellipsoid_cacheAlloc_ref,
	.cacheReset = Ellipsoid_cacheReset_ref,
	.cacheUpdateCoord = Ellipsoid_cacheUpdateCoord_ref,
	.shallowCutOracle = Ellipsoid_shallowCutOracle_ref,
	.transform = Ellipsoid_transform_ref,
        .boundingSphere = Ellipsoid_bounding_ref,
	.normal = Ellipsoid_normal,
};

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

void Ellipsoid_free(const void* o) {
   Ellipsoid* e = (Ellipsoid*)o;
   free(e->A);// includes a (and T)
   free(e);
}

void* Ellipsoid_clone(const void* o) {
   Ellipsoid* old = (Ellipsoid*)o;
   const int n = old->n;
   Ellipsoid* e = Ellipsoid_new(n);
 
   for(int i=0; i<n; i++) {
      for(int j=0; j<n; j++) {
	 e->A[i*e->line + j] = old->A[i*e->line + j];
      }
      e->a[i] = old->a[i];
   }
   
   return e;
}

void* Ellipsoid_clone_with_T(const void* o) {
   Ellipsoid* old = (Ellipsoid*)o;
   assert(old->T && "T must be allocated");
   const int n = old->n;
   Ellipsoid* e = Ellipsoid_new_with_T(n);
 
   for(int i=0; i<n; i++) {
      for(int j=0; j<n; j++) {
	      e->A[i*e->line + j] = old->A[i*e->line + j];
	      e->T[i*e->line + j] = old->T[i*e->line + j];
      }
      e->a[i] = old->a[i];
   }
   
   return e;
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
   const FT discr = b*b - 4.0*a*c;
   assert(discr >= 0);
   const FT sqrtDiscr = sqrt(discr);
   const FT aInv = 0.5/a;

   *t0 = (-b - sqrtDiscr) * aInv;
   *t1 = (-b + sqrtDiscr) * aInv;
}

void Ellipsoid_intersectCoord_ref(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache) {
   Ellipsoid* e = (Ellipsoid*)o;
   const int n = e->n;
   //FT* Az_c = (FT*)cache;

   FT* Ad = Ellipsoid_get_Ai(e,d);
   FT a = Ad[d];
   FT b = 0;
   FT c = -1.0;
   
   // do multiplications same as in eval.
   for(int i=0;i<n;i++) {
      const FT* Ai = Ellipsoid_get_Ai(e,i);
      FT Az = 0;
      for(int j=0; j<n; j++) {
         Az += Ai[j] * (x[j] - e->a[j]);
      }
      b += (i==d) * Az; // selection
      c += (x[i] - e->a[i]) * Az;
   }
   b *= 2.0;

   // find t:
   const FT discr = b*b - 4.0*a*c;
   assert(discr >= 0);
   const FT sqrtDiscr = sqrt(discr);
   const FT aInv = 0.5/a;

   *t0 = (-b - sqrtDiscr) * aInv;
   *t1 = (-b + sqrtDiscr) * aInv;
}

void Ellipsoid_intersectCoord_cached_ref(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache) {
   Ellipsoid* e = (Ellipsoid*)o;
   const int n = e->n;
   FT* Az_c = (FT*)cache;

   FT a = Ellipsoid_get_a(e,d,d);
   FT b = 2.0*Az_c[d];
   FT c = Az_c[n];
   
   // find t:
   const FT discr = b*b - 4.0*a*c;
   assert(discr >= 0);
   const FT sqrtDiscr = sqrt(discr);
   const FT aInv = 0.5/a;

   *t0 = (-b - sqrtDiscr) * aInv;
   *t1 = (-b + sqrtDiscr) * aInv;
}

int  Ellipsoid_cacheAlloc_ref(const void* o) {
   const Ellipsoid* e = (Ellipsoid*)o;
   return (e->n + 1) * sizeof(FT);
}

void Ellipsoid_cacheReset_ref(const void* o, const FT* x, void* cache) {
   const Ellipsoid* e = (Ellipsoid*)o;
   FT* c = (FT*)cache;
   const int n = e->n;
   c[n] = -1.0;
   for(int i=0;i<n;i++) {
      const FT* Ai = Ellipsoid_get_Ai(e,i);
      c[i] = 0;
      for(int j=0; j<n; j++) {
         c[i] += Ai[j] * (x[j] - e->a[j]);
      }
      c[n] += (x[i] - e->a[i]) * c[i];
   }
}

void Ellipsoid_cacheUpdateCoord_ref(const void* o, const int d, const FT dx, void* cache) {
   Ellipsoid* e = (Ellipsoid*)o;
   const int n = e->n;
   FT* c = (FT*)cache;
   c[n] += dx*(2*c[d] + dx*Ellipsoid_get_a(e,d,d));
   for(int i=0; i<n; i++) {
      c[i] += dx * Ellipsoid_get_a(e,i,d);
   }
}

bool Ellipsoid_shallowCutOracle_ref(const void* o, const Ellipsoid* e, FT* v, FT* c) {
   const Ellipsoid* this = (const Ellipsoid*)o;
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
   for(int i=0;i<n;i++) {x0[i]=this->a[i] - x0[i];}
   FT beta2 = 1.0 / (4*n*n);
   Ellipsoid_minimize(this,1, e, x0);
   Ellipsoid_minimize(this,1, e, x1);
   
   FT eval0 = Ellipsoid_eval(e,x0); 
   FT eval1 = Ellipsoid_eval(e,x1);
   
   //printf("eval: %f %f vs %f\n",eval0,eval1,beta2);
   
   //for(int i=0;i<n;i++) {printf("%.12f ",x0[i]);} printf(" - x0\n");
   //for(int i=0;i<n;i++) {printf("%.12f ",x1[i]);} printf(" - x1\n");
   if(eval0 > beta2 && eval1 > beta2) {
      if(volumeVerbose>=2) { printf("no cut found!\n"); }
      return false; // both local minima too far out
   }
   // choose better x:
   FT* x = x0;
   if(eval1 < eval0) {x = x1;}
   
   Ellipsoid_normal(this, x, v);
   *c = dotProduct(v,x, n);
   return true;
}

void Ellipsoid_transform_ref(const void* o_in, void* o_out, const Matrix* L, const FT* a, const FT beta) {
   Ellipsoid* e_in = (Ellipsoid*)o_in;
   Ellipsoid* e_out = (Ellipsoid*)o_out;
   const int n = e_in->n;
   
   // computation according to preprocess_ref
   // B' = LT * B * L
   // B'' = B' * beta^2
   
   // the matrix multiplications below are terrible, but they work 
   FT* LtB = (FT*)(aligned_alloc(32, n*n*sizeof(FT))); // align this to 32
   
   for(int i=0;i<n;i++) {
      for(int j=0;j<n;j++) {
         FT sum = 0;
	 for(int k=0;k<n;k++) {
	    FT* Lk = Matrix_get_row(L,k);
	    FT* Bk = Ellipsoid_get_Ai(e_in,k);
	    sum += Lk[i] * Bk[j];
	 }
	 LtB[i*n + j] = sum;
      }
   }
   for(int i=0;i<n;i++) {
      FT* Bi = Ellipsoid_get_Ai(e_out,i);
      for(int j=0;j<n;j++) {
         FT sum = 0;
	 for(int k=0;k<n;k++) {
	    FT* Lk = Matrix_get_row(L,k);
	    sum += LtB[i*n + k] * Lk[j];
	 }
	 Bi[j] = sum * beta * beta; // new ellipse
      }
   }

   free(LtB);
   
   // b'' = b' = L.inverse() * (a-b) / beta
   FT* ab = (FT*)(aligned_alloc(32, n*sizeof(FT))); // align this to 32
   for(int i=0;i<n;i++) {ab[i] = e_in->a[i]-a[i];}
   Matrix_L_solve(L, e_out->a, ab);
   for(int i=0;i<n;i++) {e_out->a[i] /= beta;}
   free(ab);
}

void Ellipsoid_bounding_ref(const void *B, FT *R2, FT *ori){

    const Ellipsoid *E = (Ellipsoid *) B;
    int n = E->n;
    
    // center is trivially given
    for(int i=0;i<n;i++) {ori[i] = E->a[i];}
    
    // maximize each of the 2n linear functions e_i, -e_i
    // we could also compute eigenvalues of E but this seems more difficult
    // note we can ignore center a for this
    // note the linear function e_i is maximized by sqrt(T_ii)
    // note we choose R2 = sum_{i} (max e_i - min e_i)^2 as for polytopes
    Matrix *T = Matrix_new(n,n);

    // once ellipsoid has matrix members this copying won't be necessary
    // but for now we do it, as we want to be 32byte aligned in as many functions as possible
    Matrix *A = Matrix_new(n,n);
    for (int i = 0; i < n; i++){
        FT *Ai = Matrix_get_row(A, i);
        FT *Ei = Ellipsoid_get_Ai(E, i);
        for (int j = 0; j < n; j++){
            Ai[j] = Ei[j];
        }
    }

    
    Matrix_invert_pdsym(A, T);
    for (int i = 0; i < n; i++){
        *R2 += 4*Matrix_get_row(T, i)[i];
    }

    Matrix_free(T); 
    Matrix_free(A); 
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

void Ellipsoid_normal(const void* o, const FT* x, FT* normal) {
   const Ellipsoid *e = (Ellipsoid *) o;
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

   //for(int i=0;i<n;i++) {printf("%.12f ",x[i]);} printf(" - x\n");
   //Ellipsoid_T.print(e);
   //Ellipsoid_T.print(f);
   Ellipsoid_project(e,eFac,x);
   //for(int i=0;i<n;i++) {printf("%.12f ",x[i]);} printf(" - x\n");
   
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
         if(volumeVerbose>0) { printf("Warning: taking too many steps (%d), abort minimization now! (dot: %.12f)\n",count,dot); }
	 break;
      }
   } while(dot > 0.00000001);
   
   if(volumeVerbose>=2) { printf("ellipsoid optimization: steps taken: %d, dot %.12f, eval: %.12f\n",count, dot, eval); }
   
   //for(int i=0; i<n; i++) {printf(" %.12f", x[i]);}// debug
   //printf("\n");

   free(nE);
   free(nF);
   free(nP);
}



void Ellipsoid_A_from_T(Ellipsoid* e) {
   const int n = e->n;
   if(volumeVerbose>=2) { printf("redoing A from T:\n"); }

   Matrix* L = Matrix_new(n,n);
   Matrix* Linvt = Matrix_new(n,n);
   FT* b = (FT*)(aligned_alloc(32, n*sizeof(FT))); // align this to 32
   int err = cholesky_ellipsoid(e,L);
   assert(err==0 && "no cholesky errors");
   
   for(int i=0;i<n;i++) {
      for(int j=0;j<n;j++) {b[j]=(i==j);}// unit vec
      FT* x = Matrix_get_row(Linvt, i);
      Matrix_L_solve(L, x, b);
   }
   //Matrix_print(Linvt);
   
   // A = T.inverse()
   // A = (LLt).inverse()
   // A = Lt.inverse() * L.inverse()
   for(int i=0;i<n;i++) {
      FT* Ai = Ellipsoid_get_Ai(e,i);
      for(int j=0;j<n;j++) {
	 FT* a = Matrix_get_row(Linvt, i);
	 FT* b = Matrix_get_row(Linvt, j);
         FT aij = Ai[j];
	 FT dot = dotProduct(a,b, n);
	 assert(abs(aij-dot) < 0.00000001);
	 Ai[j] = dot;
      }
   }

   free(b);
   Matrix_free(Linvt);
   Matrix_free(L);
}
