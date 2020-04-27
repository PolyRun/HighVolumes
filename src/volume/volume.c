#include "volume.h"



// dotProduct:
dotProduct_f_t dotProduct = dotProduct_ref;
#include "linalg/dotProduct/dotProduct.c"

// TODO move
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
	.transform = Ellipsoid_transform_ref,
        .boundingSphere = Ellipsoid_bounding_ref
};


preprocess_f_t preprocess_generic = preprocess_ref;

void preprocess_ref(const int n, const int bcount, const void** body_in, void** body_out, const Body_T** type, FT *det) {
   // 1. init_ellipsoid:
   //     idea: origin at 0, radius determined by min of all bodies
   printf("init_ellipsoid\n");

   // compute an enclosing ball for each body and merge them
   // TODO: is it ok if center of enclosing ball is not included in intersection of all bodies? i assume yes, the proof of theorem 3.3.9 in shallow beta-cut papersuggests that we only need that the starting ellipsoid contains all bodies

   FT R2, R2_new;
   FT *ori = (FT *) malloc(n*sizeof(FT));
   FT *ori_new = (FT *) malloc(n*sizeof(FT));
   FT *dir = (FT *) malloc(n*sizeof(FT));

   type[0]->boundingSphere(body_in[0], &R2, &ori);
   
   for (int i = 1; i < bcount; i++){
       type[i]->boundingSphere(body_in[i], &R2_new, &ori_new);

       // vector between the two origins
       for (int j = 0; j < n; j++) {
           dir[j] = ori[j] - ori_new[j];
       }
       FT dist = sqrt(dotProduct(dir, dir, n));
       FT R = sqrt(R2);
       FT R_new = sqrt(R_new);

       assert((dist < R + R_new) && "balls must intersect");
       
       // one ball contained in the other
       if (dist + R_new <= R || dist + R <= R_new) {
           // new ball contained in old ball
           if (dist + R_new <= R) {
               R2 = R2_new;
               for (int j = 0; j < n; j++) ori[j] = ori_new[j];
           }
       }
       // balls intersect in two points
       else {
           // find new center
           // this is a point on ori - ori_new at distance d1 from ori_new and d2 from ori
           // s.t. d1 + d2 = dist
           //      R2_new - d1^2 = R2 - d2^2

           FT d1 = (R2_new - R2 + dist*dist)/(2*dist);
           for (int j = 0; j < n; j++){
               ori[j] = ori_new[j] + (d1/dist) * dir[j];
           }

           // find new radius squared
           R2 = R2_new - d1*d1;
       }

   }

   printf("R2: %f\nOri: ", R2);
   for (int i = 0; i < n; i++){
       printf("%f ", ori[i]);
   }
   printf("\n");
   
   //FT R2 = 1e4; // TODO - ask sub bodies for radius, take min
   // Note: tests run reliably with 1e3
   // for 1e4 we need the periodic inverse recalculation
   // but beyond 1e5, the tests will fail
   // doing inverse adjustment always also fails the tests unfortunately...
   // I hope this problem goes away with the bounding box oracle,
   // maybe then we don't even need the inverse adjustment

   Ellipsoid* e = Ellipsoid_new_with_T(n); // origin zero
   for(int i=0; i<n; i++) {
      FT* Ai = Ellipsoid_get_Ai(e,i);
      FT* Ti = Ellipsoid_get_Ti(e,i);
      Ai[i] = 1.0/R2; // sphere with 
      Ti[i] = R2; // sphere with
      e->a[i] = ori[i];
   }
   //free(ori);
   //free(ori_new);
   //free(dir);
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
      // A' = A + fac2 * A  * Tv * Tvt * A / (1 - fac2 * Tvt * A * Tv)
      
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
      FT div = 1.0 / (1.0 - fac2*TvtATv);

      for(int i=0;i<n;i++) {
         FT* Ai = Ellipsoid_get_Ai(e,i);
         for(int j=0;j<n;j++) {
	    Ai[j] = (Ai[j] + fac2 * ATv[i]*ATv[j] * div) / zs;
	 }
      }
      
      // testing: do inverse recalculation periodically!
      if(step % 100 == 0) { // could finetune this!
         Ellipsoid_A_from_T(e);
      }

      // debug: test inverse:
      //Ellipsoid_T.print(e);
      //printf("\n");
      for(int i=0;i<n;i++){
         FT* Ai = Ellipsoid_get_Ai(e,i);
         for(int j=0;j<n;j++){
            FT sum = 0;
            for(int k=0;k<n;k++){
               FT* Tk = Ellipsoid_get_Ti(e,k);
               sum += Ai[k] * Tk[j];
            }
            //printf(" %.12f",sum);
            assert(abs(sum - 1.0*(i==j) < 0.0001));
         }
         //printf("\n");
      }
      //printf("\n");
      //assert(false);
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
   // B'' = B' * beta^2 
   // b'' = b'

   printf("Transform\n");
  
   Matrix* L = Matrix_new(n,n);
   int err = cholesky_ellipsoid(e,L);
   if (err > 0){
      printf("The input polytope is degenerate or non-existant and the volume is 0.\n");
      exit(1);		
   }
   Matrix_print(L);
   
   //printf("\nL:\n");
   //for(int i=0;i<n;i++) {
   //   FT* Li = Matrix_get_row(L,i);
   //   for(int j=0;j<n;j++) {
   //      printf(" %.12f",Li[j]);
   //   }
   //   printf("\n");
   //}
   
   // transform bodies:
    for(int b=0; b<bcount; b++) {
      type[b]->transform(body_in[b], body_out[b], L, e->a, beta);
   }

   // calculate determinant:
   //    diagonal of L matrix
   //    scaled by beta
   *det = 1;
   for (int i = 0; i < n; i++){
      *det *= Matrix_get(L,i,i);
   }
   *det /= pow(beta_r, n);

   Matrix_free(L);
   //assert(false && "fixme T");
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

FT* volume_x_ptr = NULL;
FT* volume_d_ptr = NULL;
void* volume_cache_ptr = NULL;

void volume_lib_init(const int max_n, const int max_b) {
   printf("volume_lib_init...\n");

   volume_x_ptr = (FT*)(aligned_alloc(32, max_n*sizeof(FT))); // align this to 32
   volume_d_ptr = (FT*)(aligned_alloc(32, max_n*sizeof(FT))); // align this to 32
   
   int cache_size = 1000*max_n*max_b*sizeof(FT);
   volume_cache_ptr = (aligned_alloc(32, cache_size*sizeof(FT))); // align this to 32
}

volume_f_t volume = volume_ref;

size_t pc_volume_l = 0;
size_t pc_volume_steps = 0;

FT volume_ref(const int n, const FT r0, const FT r1, const int bcount, const void** body, const Body_T** type) {
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
   FT* x = volume_x_ptr;// sample point x
   assert(x);
   for(int j=0;j<n;j++) {x[j]=0.0;}// origin

   FT* d = volume_d_ptr; // vector for random direction
   assert(d);

   // set up cache:
   int cache_size = 0;
   void* cache[bcount];
   //for(int c=0;c<bcount;c++) {
   //   cache_size += type[c]->cacheAlloc(body[c]);
   //}
   void* cache_base = volume_cache_ptr; // align this to 32
   cache_size = 0;
   for(int c=0;c<bcount;c++) {
      cache[c] = cache_base + cache_size;
      cache_size += ceil_cache(type[c]->cacheAlloc(body[c]), 1);
      // round up for allignment
      type[c]->cacheReset(body[c],x,cache[c]);
   }

   
   const int l = ceil(n*log(r1/r0) / log(2.0));
   pc_volume_l = l; // performance_counter
   pc_volume_steps = 0; // performance_counter
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
         pc_volume_steps++; // performance_counter
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

   return volume;
}



FT xyz_f1(const Polytope* body, const FT r, const int n) {return body->n;}
FT xyz_f2(const Polytope* body, const FT r, const int n) {return n;}
xyz_f_t xyz_f = xyz_f1;
