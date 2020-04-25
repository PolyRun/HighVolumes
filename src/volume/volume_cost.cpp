#include "volume_helper.hpp"
#include "volume_cost.hpp"
#include <cassert>
void xyz_f1_cost(const int n) {
   pc_stack().log(n,2*n, "bogus");
}

void dotProduct_cost_ref(const int n) {
   pc_stack().log(2*n,2*n*sizeof(FT), "dotProduct");
}
void vectorNorm_cost_ref(const int n) {
   pc_stack().log(2*n,n*sizeof(FT), "vectorNorm");
}
void Ball_intersectCoord_cost_ref(const int n) {
   {// frame for vectorNorm
      PC_Frame<vectorNorm_cost_f> frame((void*)vectorNorm);
      frame.costf()(n);
   }

   // read 1
   // div 1
   // add 4
   // mul 9
   // sqrt 1
   pc_stack().log(15,sizeof(FT), "quad. eq.");
}
void Ball_intersect_cost_ref(const int n) {
   {// frame for vectorNorm
      PC_Frame<vectorNorm_cost_f> frame((void*)vectorNorm,2); // 2 vectorNorms
      frame.costf()(n);
   }
   {// frame for dotProduct
      PC_Frame<dotProduct_cost_f> frame((void*)dotProduct);
      frame.costf()(n);
   }
   
   // div 1
   // add 4
   // mul 9
   // sqrt 1
   pc_stack().log(15,0, "quad. eq.");
}

void Polytope_intersect_cost_ref(const void* o) {
   const Polytope* p = (Polytope*)o;
   const int n = p->n;
   const int m = p->m;

   {// frame for dotProduct: m times d*ai, m times ai*x
      PC_Frame<dotProduct_cost_f> frame((void*)dotProduct, 2*m);
      frame.costf()(n);
   }

   // read m (all of b)
   // compares ??? do we count these?
   // add m
   // div m
   pc_stack().log(2*m,m*sizeof(FT), "intersect");
}
void Polytope_intersectCoord_cost_ref(const void* o) {
   const Polytope* p = (Polytope*)o;
   const int n = p->n;
   const int m = p->m;

   {// frame for dotProduct: m times ai*x
      PC_Frame<dotProduct_cost_f> frame((void*)dotProduct, m);
      frame.costf()(n);
   }

   // read m + m (all of b, ai[d])
   // compares ??? do we count these?
   // add m
   // div m
   pc_stack().log(2*m,2*m*sizeof(FT), "intersect");
}
void Polytope_intersectCoord_cached_cost_ref(const void* o) {
   const Polytope* p = (Polytope*)o;
   const int n = p->n;
   const int m = p->m;
   // read 3*m (ai[d], b, cache)
   // compares ??? do we count these?
   // add 1 * m
   // div 1 * m
   pc_stack().log(0,0, "discuss compares");
   pc_stack().log(2*m,3*m*sizeof(FT), "read cache, calculate");
}
void Polytope_cacheUpdateCoord_cost_ref(const void* o) {
   const Polytope* p = (Polytope*)o;
   const int n = p->n;
   const int m = p->m;
   // read m
   // write m
   // mul m
   // add m
   pc_stack().log(2*m,2*m*sizeof(FT), "update cached dotProduct");
}
void Polytope_cacheReset_cost_ref(const void* o) {
   const Polytope* p = (Polytope*)o;
   const int n = p->n;
   const int m = p->m;

   {// frame for dotProduct: m times ai*x
      PC_Frame<dotProduct_cost_f> frame((void*)dotProduct, m);
      frame.costf()(n);
   }

   // write m  (c)
   pc_stack().log(0,m, "write results");
}

void PolytopeT_intersect_cost_ref(const void* o) {
   const PolytopeT* p = (PolytopeT*)o;
   const int n = p->n;
   const int m = p->m;

   {// frame for dotProduct: m times d*ai, m times ai*x
      PC_Frame<dotProduct_cost_f> frame((void*)dotProduct, 2*m);
      frame.costf()(n);
   }

   // read m (all of b)
   // compares ??? do we count these?
   // add m
   // div m
   pc_stack().log(2*m,m*sizeof(FT), "intersect");
}
void PolytopeT_intersectCoord_cost_ref(const void* o) {
   const PolytopeT* p = (PolytopeT*)o;
   const int n = p->n;
   const int m = p->m;

   {// frame for dotProduct: m times ai*x
      PC_Frame<dotProduct_cost_f> frame((void*)dotProduct, m);
      frame.costf()(n);
   }

   // read m + m (all of b, ai[d])
   // compares ??? do we count these?
   // add m
   // div m
   pc_stack().log(2*m,2*m*sizeof(FT), "intersect");
}
void PolytopeT_intersectCoord_cached_cost_ref(const void* o) {
   const PolytopeT* p = (PolytopeT*)o;
   const int n = p->n;
   const int m = p->m;
   // read 3*m (A, b, cache)
   // compares ??? do we count these?
   // add 1 * m
   // div 1 * m
   pc_stack().log(0,0, "discuss compares");
   pc_stack().log(2*m,3*m*sizeof(FT), "read cache, calculate");
}
void PolytopeT_cacheUpdateCoord_cost_ref(const void* o) {
   const PolytopeT* p = (PolytopeT*)o;
   const int n = p->n;
   const int m = p->m;
   // read m
   // write m
   // mul m
   // add m
   pc_stack().log(2*m,2*m*sizeof(FT), "update cached dotProduct");
}
void PolytopeT_cacheReset_cost_ref(const void* o) {
   const PolytopeT* p = (PolytopeT*)o;
   const int n = p->n;
   const int m = p->m;
   // read n*m + n   (A, x)
   // write m  (c)
   // mul n*m
   // add n*m
   pc_stack().log(2*m*n,n*m + n + m, "recompute dotproduct");
}

void Ellipsoid_intersect_cost_ref(const void* o) {
   const Ellipsoid* e = (Ellipsoid*)o;
   const int n = e->n;
   
   // read: n^2 + n + n (all of A, all of a, all of x)
   size_t add = n + 2*n*n + 3*n + 3; // z = x-a, Az, Ad, dt*Ad, dt*Az, zt*Az, eq
   size_t mul = 2*n*n + 3*n + 6; // Az, Ad, dt*Ad, dt*Az, zt*Az, eq
   // sqrt 1
   // div 1
   pc_stack().log(add + mul + 1 + 1, (n*n + 2*n)*sizeof(FT), "read cache, calculate");
}
void Ellipsoid_intersectCoord_cost_ref(const void* o) {
   const Ellipsoid* e = (Ellipsoid*)o;
   const int n = e->n;
   
   // read: n^2 + n + n (all of A, all of a, all of x)
   size_t add = n + n*n + n + 3; // z = x-a, Az,zt*Az, eq
   size_t mul =  n*n + n + 6; // Az, zt*Az, eq
   // sqrt 1
   // div 1

   pc_stack().log(add + mul + 1 + 1, (n*n + 2*n)*sizeof(FT), "read cache, calculate");
}
void Ellipsoid_cacheUpdateCoord_cost_ref(const void* o) {
   const Ellipsoid* p = (Ellipsoid*)o;
   pc_stack().log(0,0, "no cache");
}
void Ellipsoid_cacheReset_cost_ref(const void* o) {
   pc_stack().log(0,0, "no cache");
}

void volume_cost_ref(const int n, const int bcount, const void** body, const Body_T** type) {
   pc_stack().log(0,n*sizeof(FT), "init_x");
   
   pc_stack().log(0,0, "cacheAlloc");
   pc_stack().log(0,0, "cacheReset - TODO");
   
   pc_stack().log(0,0, "Ball_volume - TODO");
   
   // number of sampling layers (steps)
   size_t l = pc_volume_l; // get from last execution
   size_t s = pc_volume_steps;

   {// frame for walk
      PC_Frame<walk_cost_f> frame((void*)walk_f,s);
      frame.costf()(n,bcount,body,type);
   }
   
   // rest of ops:
   // div 1 * l
   // mul 1 * l
   // write n * l  (reset x)
   pc_stack().log(2*l,n*l*sizeof(FT), "end of layer");

   for(int c=0;c<bcount;c++) {// body intersect
         PC_Frame<intersect_cost_f> frame((void*)type[c]->cacheReset, l);// per layer
         frame.costf()(body[c]);
   }
}
void walk_cost_ref(const int n, int bcount, const void** body, const Body_T** type) {
   int ws = walk_size;
   {// frame for walk_size loop
      PC_Frame_Base loop("loop",ws);

      pc_stack().log(0,n*sizeof(FT), "n random doubles - TODO");
      
      {// frame for Ball_intersect
         PC_Frame<Ball_intersect_cost_f> frame((void*)Ball_intersect);
         frame.costf()(n);
      }
      
      for(int c=0;c<bcount;c++) {// body intersect
         PC_Frame<intersect_cost_f> frame((void*)type[c]->intersect);
         frame.costf()(body[c]);
      }

      pc_stack().log(0,0, "random double - TODO");
      
      pc_stack().log(2*n, 3*n*sizeof(FT)," x += d*t");
   }
}

void walkCoord_cost_ref(const int n, int bcount, const void** body, const Body_T** type) {
   int ws = walk_size;
   {// frame for walk_size loop
      PC_Frame_Base loop("loop",ws);

      pc_stack().log(0,0, "random int - TODO");
      
      {// frame for Ball_intersectCoord
         PC_Frame<Ball_intersectCoord_cost_f> frame((void*)Ball_intersectCoord);
         frame.costf()(n);
      }
      
      for(int c=0;c<bcount;c++) {// body intersectCoord
         PC_Frame<intersectCoord_cost_f> frame((void*)type[c]->intersectCoord);
         frame.costf()(body[c]);
      }

      pc_stack().log(0,0, "random double - TODO");
      
      for(int c=0;c<bcount;c++) {// body intersectCoord
         PC_Frame<cacheUpdateCoord_cost_f> frame((void*)type[c]->cacheUpdateCoord);
         frame.costf()(body[c]);
      }
   }
}



