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
   // 2*m compares with +-FT_EPS and m ifs
   // add m
   // div m
   // m compares with 0.0 and m ifs
   // m compares, either t00>t or t11<t
   pc_stack().log(8*m, m*sizeof(FT), "intersect");
}

void Polytope_intersectCoord_cost_ref(const void* o) {
   const Polytope* p = (Polytope*)o;
   const int n = p->n;
   const int m = p->m;

   {// frame for dotProduct: m times ai*x
      PC_Frame<dotProduct_cost_f> frame((void*)dotProduct, m);
      frame.costf()(n);
   }

   // read 3*m (all of b, ai[d] and Aix[i])
   // 2*m compares with +-FT_EPS and m ifs
   // add m
   // div m
   // m compares with 0.0 and m ifs
   // m compares, either t00>t or t11<t
   ///////////////////////////////////////////////////////////////////////////
   // assert reads ???
   pc_stack().log(8*m, 3*m*sizeof(FT), "intersect");
}

void Polytope_intersectCoord_cached_cost_ref(const void* o) {
   const Polytope* p = (Polytope*)o;
   const int n = p->n;
   const int m = p->m;

   // read 3*m (ai[d], b, cache)
   // 2*m compares with +-FT_EPS and m ifs
   // add 1 * m
   // div 1 * m
   // m compares with 0.0 and m ifs
   // m compares, either t00>t or t11<t
   pc_stack().log(8*m, 3*m*sizeof(FT), "read cache, calculate");
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

   ////////////////////////////////////////////////////////////////////
   // read m doubles from Ai ???
   // write m  (c)
   pc_stack().log(0, 2*m, "write results");
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
   // 2*m compares with +-FT_EPS and m ifs
   // add m
   // div m
   // m compares with 0.0 and m ifs
   // m compares, either t00>t or t11<t
   pc_stack().log(8*m, m*sizeof(FT), "intersect");
}
void PolytopeT_intersectCoord_cost_ref(const void* o) {
   const PolytopeT* p = (PolytopeT*)o;
   const int n = p->n;
   const int m = p->m;
   
   pc_stack().log(0,0, "Note: early 'continue' can speed up things!");
   pc_stack().log(m*n*2,m*n*2*sizeof(FT), "dotProduct implemented locally because column-format");

   // read m + m (all of b, ai[d])
   // 2*m compares with +-FT_EPS and m ifs
   // add m
   // div m
   // m compares with 0.0 and m ifs
   // m compares, either t00>t or t11<t
   pc_stack().log(8*m, 2*m*sizeof(FT), "intersect");
}
void PolytopeT_intersectCoord_cached_cost_ref(const void* o) {
   const PolytopeT* p = (PolytopeT*)o;
   const int n = p->n;
   const int m = p->m;

   // read 3*m (A, b, cache)
   // 2*m compares with +-FT_EPS and m ifs
   // add 1 * m
   // div 1 * m
   // m compares with 0.0 and m ifs
   // m compares, either t00>t or t11<t
   pc_stack().log(8*m, 3*m*sizeof(FT), "read cache, calculate");
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
void PolytopeT_intersectCoord_cached_b_cost_ref(const void* o) {
   const PolytopeT* p = (PolytopeT*)o;
   const int n = p->n;
   const int m = p->m;

   // read 3*m (A, b, cache)
   // 2*m compares with +-FT_EPS and m ifs
   // add 1 * m
   // div 1 * m
   // m compares with 0.0 and m ifs
   // m compares, either t00>t or t11<t
   pc_stack().log(8*m, 3*m*sizeof(FT), "read cache, calculate");
   pc_stack().log(0,0, "TODO - update after impl!");
}
void PolytopeT_cacheUpdateCoord_b_cost_ref(const void* o) {
   const PolytopeT* p = (PolytopeT*)o;
   const int n = p->n;
   const int m = p->m;
   // read 2*m
   // write m
   // mul m
   // add m
   pc_stack().log(2*m,3*m*sizeof(FT), "update cached dotProduct");
   pc_stack().log(0,0, "TODO - update after impl!");
}
void PolytopeT_cacheReset_b_cost_ref(const void* o) {
   const PolytopeT* p = (PolytopeT*)o;
   const int n = p->n;
   const int m = p->m;
   // read n*m + n   (A, x)
   // write m  (c)
   // mul n*m
   // add n*m
   pc_stack().log(2*m*n,n*m + n + m, "recompute dotproduct");
   pc_stack().log(0,0, "TODO - update after impl!");
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
   size_t comp = n; // selection
   // sqrt 1
   // div 1

   pc_stack().log(add + mul + comp + 1 + 1, (n*n + 2*n)*sizeof(FT), "1 MVM, 1 VVM");
}

void Ellipsoid_intersectCoord_cached_cost_ref(const void* o) {
   const Ellipsoid* e = (Ellipsoid*)o;
   const int n = e->n;
   
   // read: 1 + 1+1 (one element of A, two elements out of cache)
   size_t add = 3; //  eq
   size_t mul = 6; // eq
   // sqrt 1
   // div 1

   pc_stack().log(add + mul + 1 + 1, 3*sizeof(FT), "read cache, calculate");
}

void Ellipsoid_intersectCoord_cached_cost_reord_fma(const void* o) {
   const Ellipsoid* e = (Ellipsoid*)o;
   const int n = e->n;
   
   // read: 1 + 1+1 (one element of A, two elements out of cache)
   // flops: see reasoning in file
   
   pc_stack().log(12, 3*sizeof(FT), "read cache, calculate");
}

void Ellipsoid_cacheUpdateCoord_cost_ref(const void* o) {
   const Ellipsoid* e = (Ellipsoid*)o;
   const int n = e->n;
   // read 2*n+1
   // write n+1
   // mul n+3
   // add n+2
   pc_stack().log(2*n+5,(3*n+2)*sizeof(FT), "update cached MVM, c");
}
void Ellipsoid_cacheReset_cost_ref(const void* o) {
   const Ellipsoid* e = (Ellipsoid*)o;
   const int n = e->n;
   // read n*n + n + n   (A, x, a)
   // write n+1  (cache)
   // mul n*n+n
   // add n*n+n
   pc_stack().log(2*n*n+2*n,(n*n + 2*n + n + 1)*sizeof(FT), "recompute MVM, c");
}


// NOTE: the actual #flops & #bytes depends on direction d
// here we return the values corresponding to the average of non-zeros per column as we expect the average to hold on average ^^
int nonzerosCSC(const PolytopeCSC *p);
int nonzerosCSC(const PolytopeCSC *p){

    // compute average non-zeros
    int nz = 0;
    for (int i = 0; i < p->col_start[p->n]; i++){
        if (p->row_idx[i] > -1){
            nz++;
        }
    }
    return nz;
}


// NOTE: the actual #flops & #bytes depends on direction d (c.f. nonzerosCSC)
void PolytopeCSC_intersect_cost_ref(const void *o){
    const PolytopeCSC *p = (PolytopeCSC *) o;
    int n = p->n;
    int m = p->m;
    int nz = nonzerosCSC(p);

    // read 4m (read and write dotd, dotx) + nz read a + m read b doubles
    // read 2n (dir, x) doubles
    // read nz ints (row_idx)
    // add 2*#non-zeros + m
    // mult 2*#non-zeros
    // div m
    // compare 3m (upper bound)
    pc_stack().log(4*nz + 5*m, (4*m + nz + m + 2*n) * sizeof(FT) + nz * sizeof(int), "intersect CSC");
}


void PolytopeCSC_mvm_cost(const PolytopeCSC *p){

    int nz = nonzerosCSC(p);

    // reads 2m (read and write result) + n (read x) + #non-zeros (read A)
    // read #non-zeros ints (for row_idx)
    // adds #non-zeros
    // mults #non-zeros
    pc_stack().log(2*nz, (2*p->m + p->n + nz) * sizeof(FT) + nz * sizeof(int), "mvm for CSC");
}

// NOTE: the actual #flops & #bytes depends on direction d (c.f. nonzerosCSC)
void PolytopeCSC_intersectCoord_cost_ref(const void *o){
    const PolytopeCSC *p = (PolytopeCSC *) o;
    int n = p->n;
    int m = p->m;
    int nz = nonzerosCSC(p);


    {// frame for mvm
        PC_Frame<mvm_cost_f> frame((void*) PolytopeCSC_mvm);
        frame.costf()(p);
    }

    // read 3*#non-zeros in col (for b, A and dotx (don't count cache here as its only for validation))
    // read #non-zeros ints (row_idx)
    // adds #non-zeros in col (1 flop each)
    // division #non-zeros in col (1 flop each)
    // comparison #non-zeros in col (1 flop each)
    pc_stack().log(3* nz/n, 3*nz/n*sizeof(FT) + nz/n *sizeof(int), "intersect Coord CSC");
    
}

// NOTE: the actual #flops & #bytes depends on direction d (c.f. nonzerosCSC)
void PolytopeCSC_intersectCoord_cached_cost_ref(const void *o){
    const PolytopeCSC *p = (PolytopeCSC *) o;
    int n = p->n;
    int m = p->m;
    int nz = nonzerosCSC(p);


    // read #non-zeros in col * (3 doubles (for Aix, b, A) + 1 int (row_idx))
    // adds non-zeros in col
    // divs non-zeros in col
    // comparisons non-zeros in col

    pc_stack().log(3*nz/n, (3 * sizeof(FT) + sizeof(int)) * nz/n, "intersectCoord_ref CSC");
}


void PolytopeCSC_intersectCoord_cached_cost_withb(const void *o){
    const PolytopeCSC *p = (PolytopeCSC *) o;
    int n = p->n;
    int m = p->m;
    int nz = nonzerosCSC(p);

    // read #non-zeros in col * (2 doubles (for A and b_Aix) + 1 int (row_idx))
    // divs #non-zeros in col
    // comparisons #non-zeros in col
    pc_stack().log(2*nz/n, (2 * sizeof(FT) + sizeof(int)) * nz/n, "intersectCoord_withb CSC");

}



// NOTE: the actual #flops & #bytes depends on direction d (c.f. nonzerosCSC)
void PolytopeCSC_cacheReset_cost_ref(const void *o){
    const PolytopeCSC *p = (PolytopeCSC *) o;
    int n = p->n;
    int m = p->m;
    int nz = nonzerosCSC(p);

    // reads 2m (read-in and write cache (assume it stays in cache between zeroing and writing)) + n (read x) + nz (read A) doubles
    // read nz ints (row_idx)
    // adds #non-zeros total
    // mults #non-zeros total
    pc_stack().log(2*nz, (2*m + n + nz) * sizeof(FT) + nz * sizeof(int), "cacheReset CSC");
}

void PolytopeCSC_cacheReset_cost_withb(const void *o){
    const PolytopeCSC *p = (PolytopeCSC *) o;
    int n = p->n;
    int m = p->m;
    int nz = nonzerosCSC(p);
    // reads 3m (read and write cache and read b) + n (read x) + nz (read A) doubles
    // read nz ints (row_idx)
    // adds #non-zeros total
    // mults #non-zeros total
    pc_stack().log(2*nz, (3*m + n + nz) * sizeof(FT) + nz * sizeof(int), "cacheReset_withb CSC");
}

// NOTE: the actual #flops & #bytes depends on direction d (c.f. nonzerosCSC)
void PolytopeCSC_cacheUpdateCoord_cost_ref(const void *o){

    const PolytopeCSC *p = (PolytopeCSC *) o;
    int n = p->n;
    int nz = nonzerosCSC(p);
    
    // read 3(read-in and write cache, read A) * #non-zeros in col * sizeof(FT) bytes
    // read #non-zeros int col ints (for row_idx)
    // mults #non-zeros in col
    // adds #non-zeros in col
    pc_stack().log(2 * nz/n, (3 * nz/n) * sizeof(FT) + nz/n * sizeof(int), "cache Update coord CSC");
    
}

// NOTE: the actual #flops & #bytes depends on direction d (c.f. nonzerosCSC)
void PolytopeCSC_cacheUpdateCoord_cost_withb(const void *o){
    PolytopeCSC_cacheUpdateCoord_cost_ref(o);    
}


void PolytopeJIT_intersect_cost_ref(const void* o) {
   const PolytopeJIT* p = (PolytopeJIT*)o;
   pc_stack().log(0,0,"TODO");
}
void PolytopeJIT_intersectCoord_cost_ref(const void* o) {
   const PolytopeJIT* p = (PolytopeJIT*)o;
   const int n = p->n;
   pc_stack().log(0,0,"init, switch case jmpq");
   
   double nzAavg = p->nzA / (double)n;
   double flops = nzAavg*2; // mul, min/max
   double data = nzAavg*1; // read cache
   pc_stack().log(flops,data*sizeof(FT),"rd cache, mul, min/max");
   pc_stack().log(0,0,"write back t0,t1");
}
void PolytopeJIT_cacheUpdateCoord_cost_ref(const void* o) {
   const PolytopeJIT* p = (PolytopeJIT*)o;
   const int n = p->n;
   pc_stack().log(0,0,"switch case jmpq");

   double nzAavg = p->nzA / (double)n;
   double flops = nzAavg*2; // fmadd
   double data = nzAavg*2; // rd/wr cache
   pc_stack().log(flops,data*sizeof(FT),"rd/wr cache, fmadd");
}
void PolytopeJIT_cacheReset_cost_ref(const void* o) {
   const PolytopeJIT* p = (PolytopeJIT*)o;
   pc_stack().log(0,0,"TODO");
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



