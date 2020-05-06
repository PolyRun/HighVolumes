#include "volume_helper.hpp"
#include "volume_cost.hpp"
#include <cassert>

void xyz_f1_cost(const int n) {
   pc_stack().log(n,2*n, "bogus");
}

void dotProduct_cost_ref(const int n) {
   pc_stack().log(2*n,2*n*sizeof(FT), "dotProduct");
}

void squaredNorm_cost_ref(const int n) {
   pc_stack().log(2*n, n*sizeof(FT), "squaredNorm");
}

void Ball_intersectCoord_cost_ref(const int n) {

   // frame for squaredNorm
   {
      PC_Frame<squaredNorm_cost_f> frame((void*) squaredNorm);
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

   // frame for squaredNorm
   {
      PC_Frame<squaredNorm_cost_f> frame((void*) squaredNorm, 2); // 2 squaredNorms
      frame.costf()(n);
   }
   {// frame for dotProduct
      PC_Frame<dotProduct_cost_f> frame((void*) dotProduct);
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
   // read 2*m
   // write m
   // mul m
   // add m
   pc_stack().log(2*m,3*m*sizeof(FT), "update cached dotProduct");
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
   
   pc_stack().log(0,0, "Note: early 'continue' can speed up things!");
   pc_stack().log(m*n*2,m*n*2*sizeof(FT), "dotProduct implemented locally because column-format");

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
   // read 2*m
   // write m
   // mul m
   // add m
   pc_stack().log(2*m,3*m*sizeof(FT), "update cached dotProduct");
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
   pc_stack().log(add + mul + 1 + 1, (n*n + 2*n)*sizeof(FT), "2 MVM (parallel), some VVM");
}
void Ellipsoid_intersectCoord_cost_ref(const void* o) {
   const Ellipsoid* e = (Ellipsoid*)o;
   const int n = e->n;
   
   // read: n^2 + n + n (all of A, all of a, all of x)
   size_t add = n + n*n + n + 3; // z = x-a, Az,zt*Az, eq
   size_t mul =  n*n + n + 6; // Az, zt*Az, eq
   // sqrt 1
   // div 1

   pc_stack().log(add + mul + 1 + 1, (n*n + 2*n)*sizeof(FT), "1 MVM, 1 VVM");
}

void Ellipsoid_intersectCoord_cached_cost_ref(const void* o) {
   const Ellipsoid* e = (Ellipsoid*)o;
   const int n = e->n;
   
   // read: n^2 + n + n (all of A, all of a, all of x, cache)
   size_t add = n + n + 3; // z = x-a, eq
   size_t mul = n + 6; // Az, eq
   // sqrt 1
   // div 1

   pc_stack().log(add + mul + 1 + 1, (n*n + 3*n)*sizeof(FT), "read cache, calculate");
}

void Ellipsoid_cacheUpdateCoord_cost_ref(const void* o) {
   const Ellipsoid* e = (Ellipsoid*)o;
   const int n = e->n;
   // read 2*n
   // write n
   // mul n
   // add n
   pc_stack().log(2*n,3*n*sizeof(FT), "update cached MVM");
}
void Ellipsoid_cacheReset_cost_ref(const void* o) {
   const Ellipsoid* e = (Ellipsoid*)o;
   const int n = e->n;
   // read n*n + n   (A, x)
   // write m  (c)
   // mul n*n
   // add n*n
   pc_stack().log(2*n*n,n*n + n + n, "recompute MVM");
}

void volume_cost_ref(const int n, const int bcount, const void** body, const Body_T** type) {

   pc_stack().log(0, n*sizeof(FT), "init_x");
   
   pc_stack().log(0, 0, "cacheAlloc");
   pc_stack().log(0, 0, "cacheReset - TODO");
   
   pc_stack().log(0, 0, "Ball_volume - TODO");
   
   // number of sampling layers (steps)
   size_t l = pc_volume_l; // get from last execution
   size_t s = pc_volume_steps;

   // frame for steps loop
   {
      PC_Frame_Base loop("steps",s);
      {
         PC_Frame<walk_cost_f> frame((void*) walk_f);
         frame.costf()(n, bcount, body, type);
      }

      // frame for squaredNorm: x
      {
         PC_Frame<squaredNorm_cost_f> frame((void*) squaredNorm); // squaredNorm
         frame.costf()(n);
      }
      
      // log 2
      // div 2
      // mul 2
      pc_stack().log(6, 0, "find layer: logs!");
   }
   
   // rest of ops:
   // div 1 * l
   // mul 1 * l
   // read n * l and write n * l  (reset x)
   pc_stack().log(2*l, 2*n*l*sizeof(FT), "end of layer");

   // body intersect
   for (int c = 0; c < bcount; c++) {
      // per layer
      PC_Frame<cacheReset_cost_f> frame((void*) type[c]->cacheReset, l);
      frame.costf()(body[c]);
   }

}

void walk_cost_ref(const int n, int bcount, const void** body, const Body_T** type) {

   int ws = walk_size;

   // frame for walk_size loop
   {
      PC_Frame_Base loop("loop",ws);

      pc_stack().log(0,n*sizeof(FT), "n random doubles - TODO");
      
      {// frame for Ball_intersect
         PC_Frame<Ball_intersect_cost_f> frame((void*) Ball_intersect);
         frame.costf()(n);
      }
      
      // body intersect
      for(int c = 0; c < bcount; c++) {
         PC_Frame<intersect_cost_f> frame((void*) type[c]->intersect);
         frame.costf()(body[c]);
      }

      pc_stack().log(0,0, "random double - TODO");
      
      pc_stack().log(2*n, 3*n*sizeof(FT)," x += d*t");
   }

}

void walkCoord_cost_ref(const int n, int bcount, const void** body, const Body_T** type) {

   int ws = walk_size;

   // frame for walk_size loop
   {
      PC_Frame_Base loop("loop", ws);

      pc_stack().log(0, 0, "random int - TODO");
      
      // frame for Ball_intersectCoord
      {
         PC_Frame<Ball_intersectCoord_cost_f> frame((void*) Ball_intersectCoord);
         frame.costf()(n);
      }
      
      // body intersectCoord
      for (int c = 0; c < bcount; c++) {
         PC_Frame<intersectCoord_cost_f> frame((void*) type[c]->intersectCoord);
         frame.costf()(body[c]);
      }

      pc_stack().log(0, 0, "random double - TODO");

      // Reading and writing x[dd] with one add in between
      pc_stack().log(1, 2, "x[dd] += t;");
      
      // body intersectCoord
      for(int c = 0; c < bcount; c++) {
         PC_Frame<cacheUpdateCoord_cost_f> frame((void*) type[c]->cacheUpdateCoord);
         frame.costf()(body[c]);
      }
   }
   
}



