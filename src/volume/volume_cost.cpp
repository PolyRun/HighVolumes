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

void Random_int_cost_ref(const void* o){
    // 1 read (status)
    // 1 write (status)
    // 3 shifts
    // 1 bit-wise and
    // 3 bit-wise xor
    pc_stack().log(7, 2*sizeof(int), "random int");
}

void Random_int_in_range_cost_ref(const void* o){

    {// frame for random int
        PC_Frame<random_int_cost_f> frame((void*) prng_get_random_int);
        frame.costf()((NULL));
    }

    // 1 mod
    // 3 add
    pc_stack().log(4, 0, "random int_in_range");
}

void Random_double_in_range_cost_ref(const void* o){

    {// frame for random int
        PC_Frame<random_int_cost_f> frame((void*) prng_get_random_int);
        frame.costf()((NULL));
    }

    // 1 div
    // 3 mul
    // 2 add
    pc_stack().log(5, 0, "random double_in_range");
}
void sr_rand256d_cost_ref(){
    {// frame for random int
        PC_Frame<random_int_cost_f> frame((void*) prng_get_random_int,4);
        frame.costf()((NULL));
    }

    // 1 and
    // 1 shift
    // 1 or
    // 1 sub
    //  all times 4
    pc_stack().log(4*4, 0, "rand256d");
}

void Random_double_0_1_cost_ref(const void* o){

    {// frame for random int
        PC_Frame<random_int_cost_f> frame((void*) prng_get_random_int);
        frame.costf()((NULL));
    }

    // 1 div
    pc_stack().log(1, 0, "random double_0_1");
}

void Random_double_normal_cost_ref(const void* o){

    {// frame for random double_0_1
        PC_Frame<random_double_0_1_cost_f> frame((void*) prng_get_random_double_0_1,2); // 2 random doubles_0_1
        frame.costf()(NULL);
    }

    // 3 functions (log, sqrt, log)
    // 4 mul

    pc_stack().log(7, 0, "random double_normal");
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

void Ball_intersectCoord_cached_cost_ref(const int n) {
   // read cache
   // read 1
   // div 1
   // add 4
   // mul 9
   // sqrt 1
   pc_stack().log(15,2*sizeof(FT), "quad. eq.");
}
void Ball_intersectCoord_cached4_cost_ref(const int n) {
   // read cache
   // read 1
   // add 4
   // mul 6
   // sqrt 1
   pc_stack().log(12*4,4*2*sizeof(FT), "quad. eq.");
}
void Ball_intersectCoord_cached8_cost_ref(const int n) {
   // read cache
   // read 1
   // div 1
   // add 4
   // mul 6
   // sqrt 1
   pc_stack().log(12*8,8*2*sizeof(FT), "quad. eq.");
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

   /*
Note: the block

      if(dai < 0.0) {
         t00 = (t00>t)?t00:t; // max
      } else {
         t11 = (t11<t)?t11:t; // min
      }

gets translated to something like: 
c3df0:	c5 fb 5d 6c 24 20    	vminsd 0x20(%rsp),%xmm0,%xmm5
c3df6:	c5 fb 11 6c 24 20    	vmovsd %xmm5,0x20(%rsp)
...
c3e7e:	c5 f9 2f e1          	vcomisd %xmm1,%xmm4
...
c3e86:	0f 86 64 ff ff ff    	jbe    c3df0 <Polytope_intersect_ref+0x60>
c3e8c:	c5 fb 5f 7c 24 28    	vmaxsd 0x28(%rsp),%xmm0,%xmm7
...
c3e96:	c5 fb 11 7c 24 28    	vmovsd %xmm7,0x28(%rsp)

i.e. 1 comparison and 1 max (min respectively) per iteration

    */
   
   // read m (all of b)
   // (1) 2*m compares with +-FT_EPS
   // (2) add m
   // (3) div m
   // (4) m compares with 0.0 and m ifs (MB: i think this should be only 1*m, c.f. assembly)
   // (5) m compares, either t00>t or t11<t (MB: i.e. either min or max)
   // MB: this gives (1)+(2)+(3)+(4)+(5) = (2+1+1+1+1) * m 
   pc_stack().log(6*m, m*sizeof(FT), "intersect");
}

void Polytope_intersectCoord_cost_ref(const void* o) {
   const Polytope* p = (Polytope*)o;
   const int n = p->n;
   const int m = p->m;

   {// frame for dotProduct: m times ai*x
      PC_Frame<dotProduct_cost_f> frame((void*)dotProduct, m);
      frame.costf()(n);
   }

   // read 2*m (all of b, ai[d])
   // 2*m compares with +-FT_EPS and m ifs
   // add m
   // div m
   // m compares with 0.0 and m ifs -> c.f. intersect_ref -> only 1*m
   // m compares, either t00>t or t11<t
   ///////////////////////////////////////////////////////////////////////////
   // assert reads ???
   pc_stack().log(6*m, 2*m*sizeof(FT), "intersect");
}

void Polytope_intersectCoord_cached_cost_ref(const void* o) {
   const Polytope* p = (Polytope*)o;
   const int n = p->n;
   const int m = p->m;

   // read 3*m (ai[d], b, cache)
   // 2*m compares with +-FT_EPS
   // add 1 * m
   // div 1 * m
   // m compares with 0.0 and m ifs -> MB: as in intersect_ref -> 1*m
   // m compares, either t00>t or t11<t
   pc_stack().log(6*m, 3*m*sizeof(FT), "read cache, calculate");
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

   // read and write m from c (read before write!)
   pc_stack().log(0, 2*m*sizeof(FT), "write results");
}

void PolytopeT_intersect_cost_ref(const void* o) {
   const PolytopeT* p = (PolytopeT*)o;
   const int n = p->n;
   const int m = p->m;

   pc_stack().log(4*m*n, (m*n+2*n)*sizeof(FT), "Fuzed dot products d*a and x*a");

   // read m (all of b)
   // 2*m compares with +-FT_EPS
   // add m
   // div m
   // m compares with 0.0
   // m compares, either t00>t or t11<t
   pc_stack().log(6*m, m*sizeof(FT), "intersect");
}
void PolytopeT_intersectCoord_cost_ref(const void* o) {
   const PolytopeT* p = (PolytopeT*)o;
   const int n = p->n;
   const int m = p->m;
   
   pc_stack().log(0,0, "Note: early 'continue' can speed up things!");
   pc_stack().log(m*n*2,(m*n+n)*sizeof(FT), "dotProduct implemented locally because column-format");

   // read m (all of b)
   // 2*m compares with +-FT_EPS
   // add m
   // div m
   // m compares with 0.0
   // m compares, either t00>t or t11<t
   pc_stack().log(6*m, m*sizeof(FT), "intersect");
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
   pc_stack().log(6*m, 3*m*sizeof(FT), "read cache, calculate");
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
   // read + write m  (c) MB: also need to read c before write
   // mul n*m
   // add n*m
   pc_stack().log(2*m*n,(n*m + n + 2*m)*sizeof(FT), "recompute dotproduct");
}
void PolytopeT_intersectCoord_cached_b_cost_ref(const void* o) {
   const PolytopeT* p = (PolytopeT*)o;
   const int n = p->n;
   const int m = p->m;

   // read 2*m (A, cache)
   // 2*m compares with +-FT_EPS
   // div 1 * m 
   // m compares with 0.0
   // m compares, either t00>t or t11<t
   pc_stack().log(5*m, 2*m*sizeof(FT), "read cache, calculate");
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
}

void PolytopeT_cacheUpdateCoord_b_cost_vec(const void* o) {
   PolytopeT_cacheUpdateCoord_b_cost_ref(o);
}

void PolytopeT_cacheReset_b_cost_ref(const void* o) {
   const PolytopeT* p = (PolytopeT*)o;
   const int n = p->n;
   const int m = p->m;
   // read n*m + n + m + 2*m   (A, x, b, rw c)
   // mul n*m
   // add n*m
   pc_stack().log(2*m*n,(n*m + n + 3*m)*sizeof(FT), "recompute dotproduct with b");
}

void PolytopeT_cacheReset_b_cost_vec(const void* o) {
   PolytopeT_cacheReset_b_cost_ref(o);
}
void PolytopeT_intersectCoord4_cost_ref(const void* o) {
   const PolytopeT* p = (PolytopeT*)o;
   const int n = p->n;
   const int m = p->m;
   
   // read m: Ainv column
   // read 4*m: cache
   // cmp: 2*m*4 (non-zero check included)
   // mul: 4*m
   // max/min: 2*m*4
   pc_stack().log(20*m, 5*m*sizeof(FT), "intersect");
}
void PolytopeT_intersectCoord8_cost_ref(const void* o) {
   const PolytopeT* p = (PolytopeT*)o;
   const int n = p->n;
   const int m = p->m;
   
   // read m: Ainv column
   // read 8*m: cache
   // cmp: 2*m*8 (non-zero check included)
   // mul: m*8
   // max/min: 2*m*8
   pc_stack().log(40*m, 9*m*sizeof(FT), "intersect");
}
void PolytopeT_cacheUpdateCoord4_cost_ref(const void* o) {
   const PolytopeT* p = (PolytopeT*)o;
   const int n = p->n;
   const int m = p->m;
   // read m: clumn of A
   // read m*4: cache
   // write m*4: cache
   // mul m*4
   // add m*4
   pc_stack().log(8*m,9*m*sizeof(FT), "update cached dotProduct");
}
void PolytopeT_cacheUpdateCoord8_cost_ref(const void* o) {
   const PolytopeT* p = (PolytopeT*)o;
   const int n = p->n;
   const int m = p->m;
   // read m: clumn of A
   // read m*8: cache
   // write m*8: cache
   // mul m*8
   // add m*8
   pc_stack().log(16*m,17*m*sizeof(FT), "update cached dotProduct");
}
void PolytopeT_cacheReset4_cost_ref(const void* o) {
   const PolytopeT* p = (PolytopeT*)o;
   const int n = p->n;
   const int m = p->m;
   // read n*m + 4*n + m + 2*m*4   (A, x, b, rw c)
   // mul n*m*4
   // add n*m*4
   pc_stack().log(8*m*n,(n*m + 4*n + 9*m)*sizeof(FT), "recompute dotproduct with b");
}
void PolytopeT_cacheReset8_cost_ref(const void* o) {
   const PolytopeT* p = (PolytopeT*)o;
   const int n = p->n;
   const int m = p->m;
   // read n*m + 8*n + m + 2*m*8   (A, x, b, rw c)
   // mul n*m*8
   // add n*m*8
   pc_stack().log(16*m*n,(n*m + 8*n + 17*m)*sizeof(FT), "recompute dotproduct with b");
}


void Ellipsoid_intersect_cost_ref(const void* o) {
   const Ellipsoid* e = (Ellipsoid*)o;
   const int n = e->n;
   
   // read: n^2 + n + n + n (all of A, all of a, all of x, all of d)
   size_t add = 4*n + 3*n*n + 3;  
   size_t mul = 2*n*n + 3*n + 6; 
   // sqrt 1
   // div 1
   pc_stack().log(add + mul + 1 + 1, (n*n + 3*n)*sizeof(FT), "2 MVM (parallel), some VVM");
}
void Ellipsoid_intersectCoord_cost_ref(const void* o) {
   const Ellipsoid* e = (Ellipsoid*)o;
   const int n = e->n;
   
   // read: n^2 + n + n (all of A, all of a, all of x)
   size_t add = 2*n + 2*n*n + 3;
   size_t mul =  n*n + n + 6; 
   // sqrt 1
   // div 1
   pc_stack().log(add + mul + 1 + 1, (n*n + 2*n)*sizeof(FT), "1 MVM, 1 VVM");
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
   // flops: see reasoning in file
   Ellipsoid_intersectCoord_cached_cost_ref(o);
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
   // rd/write 2*n+2  (cache)
   // mul n*n+n
   // add 2*n*n+2*n
   pc_stack().log(3*n*n+3*n,(n*n + 4*n + 2)*sizeof(FT), "recompute MVM, c");
}




// NOTE: the actual #flops & #bytes depends on direction d (c.f. nonzerosCSC)
void PolytopeCSC_intersect_cost_ref(const void *o){
    const PolytopeCSC *p = (PolytopeCSC *) o;
    int n = p->n;
    int m = p->m;
    int nz = nonzerosCSC(p);

    // nz read A + m read b doubles
    // read 2n (dir, x) doubles
    // read nz ints (row_idx)
    // add 2*#non-zeros + m
    // mult 2*#non-zeros
    // div m
    // compare 3m (upper bound)
    // max/min m
    pc_stack().log(4*nz + 6*m, (nz + m + 2*n) * sizeof(FT) + nz * sizeof(int), "intersect CSC");
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

    // read 2*#non-zeros in col (for b, A (don't count cache here as its only for validation))
    // read #non-zeros ints (row_idx)
    // adds #non-zeros in col (1 flop each)
    // division #non-zeros in col (1 flop each)
    // comparison #non-zeros in col (1 flop each)
    // min/max 1 flop each
    pc_stack().log(4* nz/n, 2*nz/n*sizeof(FT) + nz/n *sizeof(int), "intersect Coord CSC");
    
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
    // min/max
    pc_stack().log(4*nz/n, (3 * sizeof(FT) + sizeof(int)) * nz/n, "intersectCoord_ref CSC");
}


void PolytopeCSC_intersectCoord_cached_cost_withb(const void *o){
    const PolytopeCSC *p = (PolytopeCSC *) o;
    int n = p->n;
    int m = p->m;
    int nz = nonzerosCSC(p);

    // read #non-zeros in col * (2 doubles (for A and b_Aix) + 1 int (row_idx))
    // divs #non-zeros in col
    // comparisons #non-zeros in col
    // min/max
    pc_stack().log(3*nz/n, (2 * sizeof(FT) + sizeof(int)) * nz/n, "intersectCoord_withb CSC");
}


void PolytopeCSC_intersectCoord_cached_cost_vec(const void *o){
    const PolytopeCSC *p = (PolytopeCSC *) o;
    int n = p->n;
    int m = p->m;
    int nz = nonzerosCSC(p);

    // read #non-zeros in col * (2 doubles (for A and b_Aix) + 1 int (row_idx))
    // divs #non-zeros in col
    // min and max each #non-zeros in col
    // comparison #non-zeros in col
    pc_stack().log(4*nz/n, (2 * sizeof(FT) + sizeof(int)) * nz/n, "intersectCoord_withb CSC");
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


void PolytopeCSC_intersectCoord4_cost_ref(const void* o) {
    const PolytopeCSC *p = (PolytopeCSC *) o;
    int n = p->n;
    int m = p->m;
    int nz = nonzerosCSC(p);

    // read #non-zeros in col * (5 doubles (for A and 4*b_Aix) + 1 int (row_idx))
    // 4*mul #non-zeros in col
    // 4*min and 4*max each #non-zeros in col
    // 4*comparison #non-zeros in col
    pc_stack().log(16*nz/n, (5 * sizeof(FT) + sizeof(int)) * nz/n, "intersectCoord4 CSC");
}
void PolytopeCSC_intersectCoord8_cost_ref(const void* o) {
    const PolytopeCSC *p = (PolytopeCSC *) o;
    int n = p->n;
    int m = p->m;
    int nz = nonzerosCSC(p);

    // read #non-zeros in col * (9 doubles (for A and 8*b_Aix) + 1 int (row_idx))
    // 8*mul #non-zeros in col
    // 8*min and 8*max each #non-zeros in col
    // 8*comparison #non-zeros in col
    pc_stack().log(32*nz/n, (9 * sizeof(FT) + sizeof(int)) * nz/n, "intersectCoord8 CSC");
}
void PolytopeCSC_cacheUpdateCoord4_cost_ref(const void* o) {
    const PolytopeCSC *p = (PolytopeCSC *) o;
    int n = p->n;
    int nz = nonzerosCSC(p);
    
    // read 8+1(read-in and write cache*4, read A) * #non-zeros in col * sizeof(FT) bytes
    // read #non-zeros int col ints (for row_idx)
    // 4*mults #non-zeros in col
    // 4*adds #non-zeros in col
    pc_stack().log(8 * nz/n, (9 * nz/n) * sizeof(FT) + nz/n * sizeof(int), "cache Update coor4d CSC");
}
void PolytopeCSC_cacheUpdateCoord8_cost_ref(const void* o) {
    const PolytopeCSC *p = (PolytopeCSC *) o;
    int n = p->n;
    int nz = nonzerosCSC(p);
    
    // read 16+1(read-in and write cache*8, read A) * #non-zeros in col * sizeof(FT) bytes
    // read #non-zeros int col ints (for row_idx)
    // 8*mults #non-zeros in col
    // 8*adds #non-zeros in col
    pc_stack().log(16 * nz/n, (17 * nz/n) * sizeof(FT) + nz/n * sizeof(int), "cache Update coord8 CSC");
}
void PolytopeCSC_cacheReset4_cost_ref(const void* o) {
    const PolytopeCSC *p = (PolytopeCSC *) o;
    int n = p->n;
    int m = p->m;
    int nz = nonzerosCSC(p);
    // reads (8+1)m (read and write cache*4 and read b) + 4*n (read x) + nz (read A) doubles
    // read nz ints (row_idx)
    // 4*adds #non-zeros total
    // 4*mults #non-zeros total
    pc_stack().log(8*nz, (9*m + 4*n + nz) * sizeof(FT) + nz * sizeof(int), "cacheReset4");
}
void PolytopeCSC_cacheReset8_cost_ref(const void* o) {
    const PolytopeCSC *p = (PolytopeCSC *) o;
    int n = p->n;
    int m = p->m;
    int nz = nonzerosCSC(p);
    // reads (16+1)m (read and write cache*8 and read b) + 8*n (read x) + nz (read A) doubles
    // read nz ints (row_idx)
    // 8*adds #non-zeros total
    // 8*mults #non-zeros total
    pc_stack().log(16*nz, (17*m + 8*n + nz) * sizeof(FT) + nz * sizeof(int), "cacheReset4");
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
   double data = nzAavg*2; // read cache + A
   pc_stack().log(flops,data*sizeof(FT),"rd cache+const, mul, min/max");
}
void PolytopeJIT_cacheUpdateCoord_cost_ref(const void* o) {
   const PolytopeJIT* p = (PolytopeJIT*)o;
   const int n = p->n;
   pc_stack().log(0,0,"switch case jmpq");

   double nzAavg = p->nzA / (double)n;
   double flops = nzAavg*2; // fmadd
   double data = nzAavg*3; // rd/wr cache
   pc_stack().log(flops,data*sizeof(FT),"rd/wr cache, rd const, fmadd");
}
void PolytopeJIT_cacheReset_cost_ref(const void* o) {
   const PolytopeJIT* p = (PolytopeJIT*)o;
   // read A, b, x, rw b
   // fmsub for nz
   pc_stack().log(p->nzA*2,(p->nzA + 3*p->m + p->n)*sizeof(FT),"recompute cache, dot prod");
}

void PolytopeJIT_intersectCoord4_cost_ref(const void* o) {
   const PolytopeJIT* p = (PolytopeJIT*)o;
   const int n = p->n;
   pc_stack().log(0,0,"init, switch case jmpq");
   
   double nzAavg = p->nzA / (double)n;
   double flops = nzAavg*2*4; // mul, min/max  *4
   double data = nzAavg*5; // read cache*4 + A
   pc_stack().log(flops,data*sizeof(FT),"rd cache+const, mul, min/max");
}
void PolytopeJIT_cacheUpdateCoord4_cost_ref(const void* o) {
   const PolytopeJIT* p = (PolytopeJIT*)o;
   const int n = p->n;
   pc_stack().log(0,0,"switch case jmpq");

   double nzAavg = p->nzA / (double)n;
   double flops = nzAavg*2*4; // fmadd
   double data = nzAavg*9; // rd/wr cache*8, A
   pc_stack().log(flops,data*sizeof(FT),"rd/wr cache, rd const, fmadd");
}
void PolytopeJIT_cacheReset4_cost_ref(const void* o) {
   const PolytopeJIT* p = (PolytopeJIT*)o;
   // read A, b, x, rw c
   // fmsub for nz
   pc_stack().log(p->nzA*2*4,(p->nzA + (1+8)*p->m + 4*p->n)*sizeof(FT),"recompute cache, dot prod");
}

void PolytopeJIT_intersectCoord8_cost_ref(const void* o) {
   const PolytopeJIT* p = (PolytopeJIT*)o;
   const int n = p->n;
   pc_stack().log(0,0,"init, switch case jmpq");
   
   double nzAavg = p->nzA / (double)n;
   double flops = nzAavg*2*8; // mul, min/max  *8
   double data = nzAavg*8; // read cache*8 + A
   pc_stack().log(flops,data*sizeof(FT),"rd cache+const, mul, min/max");
}
void PolytopeJIT_cacheUpdateCoord8_cost_ref(const void* o) {
   const PolytopeJIT* p = (PolytopeJIT*)o;
   const int n = p->n;
   pc_stack().log(0,0,"switch case jmpq");

   double nzAavg = p->nzA / (double)n;
   double flops = nzAavg*2*8; // fmadd
   double data = nzAavg*17; // rd/wr cache*16, A
   pc_stack().log(flops,data*sizeof(FT),"rd/wr cache, rd const, fmadd");
}
void PolytopeJIT_cacheReset8_cost_ref(const void* o) {
   const PolytopeJIT* p = (PolytopeJIT*)o;
   // read A, b, x, rw c
   // fmsub for nz
   pc_stack().log(p->nzA*2*8,(p->nzA + (1+16)*p->m + 8*p->n)*sizeof(FT),"recompute cache, dot prod");
}

void volume_cost_ref(const int n, const int bcount, const void** body, const Body_T** type) {

   pc_stack().log(0, 2*n*sizeof(FT), "init_x");
   
   pc_stack().log(0, 0, "cacheAlloc");
   
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
      // ceil
      // comp
      pc_stack().log(8, 0, "find layer: logs!");
   }
   
   // rest of ops:
   // div 1 * l
   // mul 1 * l
   // read n * l and write n * l  (reset x)
   // mul n * l
   pc_stack().log(2*l + n*l, 2*n*l*sizeof(FT), "end of layer");

   // body intersect
   for (int c = 0; c < bcount; c++) {
      // per layer
      PC_Frame<cacheReset_cost_f> frame((void*) type[c]->cacheReset, l+1);// +1 for init
      frame.costf()(body[c]);
   }

}

void volume_coord_1_cost_ref(const int n, const int bcount, const void** body, const Body_T** type) {

   pc_stack().log(0, 2*n*sizeof(FT), "init_x");
   
   pc_stack().log(0, 0, "cacheAlloc");
   
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

      pc_stack().log(0, sizeof(FT), "read out cached squaredNorm");
      
      // log 2
      // div 2
      // mul 2
      // ceil
      // comp
      pc_stack().log(8, 0, "find layer: logs!");
   }
   
   // rest of ops:
   // div 1 * l
   // mul 1 * l
   // read n * l and write n * l  (reset x)
   // mul n * l
   pc_stack().log(2*l + n*l, 2*n*l*sizeof(FT), "end of layer");

   // body intersect
   for (int c = 0; c < bcount; c++) {
      // per layer
      PC_Frame<cacheReset_cost_f> frame((void*) type[c]->cacheReset, l+1);// +1 for init
      frame.costf()(body[c]);
   }

}

void volume_coord_4_cost_ref(const int n, const int bcount, const void** body, const Body_T** type) {
   pc_stack().log(0, 4*2*n*sizeof(FT), "init_x");
   
   pc_stack().log(0, 0, "cacheAlloc");
   
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

      pc_stack().log(0, 4*sizeof(FT), "read out cached squaredNorm");
      
      // log 2
      // div 2
      // mul 2
      // ceil
      // comp
      // * 4
      pc_stack().log(4*8, 0, "find layer: logs!");
   }
   
   // rest of ops:
   // div 1 * l
   // mul 1 * l
   //
   // read n * l and write n * l  (reset x)     *4
   // mul n * l   *4
   pc_stack().log(2*l + n*l*4, 4*2*n*l*sizeof(FT), "end of layer");

   // body intersect
   for (int c = 0; c < bcount; c++) {
      // per layer
      PC_Frame<cacheReset_cost_f> frame((void*) type[c]->cacheReset4, l+1);// +1 for init
      frame.costf()(body[c]);
   }


}
void volume_coord_8_cost_ref(const int n, const int bcount, const void** body, const Body_T** type) {
   pc_stack().log(0, 8*2*n*sizeof(FT), "init_x");
   
   pc_stack().log(0, 0, "cacheAlloc");
   
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

      pc_stack().log(0, 8*sizeof(FT), "read out cached squaredNorm");
      
      // log 2
      // div 2
      // mul 2
      // ceil
      // comp
      // * 8
      pc_stack().log(8*8, 0, "find layer: logs!");
   }
   
   // rest of ops:
   // div 1 * l
   // mul 1 * l
   //
   // read n * l and write n * l  (reset x)     *4
   // mul n * l   *4
   pc_stack().log(2*l + n*l*8, 8*2*n*l*sizeof(FT), "end of layer");

   // body intersect
   for (int c = 0; c < bcount; c++) {
      // per layer
      PC_Frame<cacheReset_cost_f> frame((void*) type[c]->cacheReset8, l+1);// +1 for init
      frame.costf()(body[c]);
   }
}

void walk_cost_ref(const int n, int bcount, const void** body, const Body_T** type) {

   int ws = walk_size;

   // frame for walk_size loop
   {
      PC_Frame_Base loop("loop",ws);

      pc_stack().log(0,2*n*sizeof(FT), "rw n random doubles");
      {// frame for random double_normal
         PC_Frame<random_double_normal_cost_f> frame((void*) prng_get_random_double_normal,n);
         frame.costf()(NULL);
      }
      
      {// frame for Ball_intersect
         PC_Frame<Ball_intersect_cost_f> frame((void*) Ball_intersect);
         frame.costf()(n);
      }
      
      // body intersect
      for(int c = 0; c < bcount; c++) {
         PC_Frame<intersect_cost_f> frame((void*) type[c]->intersect);
         frame.costf()(body[c]);
      }
      pc_stack().log(2*bcount,0,"Update min/max intersection point for last body.");

      {// frame for random double_in_range
         PC_Frame<random_double_in_range_cost_f> frame((void*) prng_get_random_double_in_range);
         frame.costf()(NULL);
      }
      
      pc_stack().log(2*n, 2*n*sizeof(FT)," x += d*t");// don't recount d here
   }

}

void walkCoord_cost_ref(const int n, int bcount, const void** body, const Body_T** type) {

   int ws = walk_size;

   // frame for walk_size loop
   {
      PC_Frame_Base loop("loop", ws);

      {// frame for random int_in_range
        PC_Frame<random_int_cost_f> frame((void*) prng_get_random_int_in_range);
        frame.costf()((NULL));
    }
      
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
      pc_stack().log(2*bcount,0,"Update min/max intersection point for last body.");

      {// frame for random double_in_range
         PC_Frame<random_double_in_range_cost_f> frame((void*) prng_get_random_double_in_range);
         frame.costf()(NULL);
      }

      // Reading and writing x[dd] with one add in between
      pc_stack().log(1, 2*sizeof(FT), "x[dd] += t;");
      
      // body intersectCoord
      for(int c = 0; c < bcount; c++) {
         PC_Frame<cacheUpdateCoord_cost_f> frame((void*) type[c]->cacheUpdateCoord);
         frame.costf()(body[c]);
      }
   }
   
}

void walkCoord_coord_1_cost_ref(const int n, int bcount, const void** body, const Body_T** type) {

   int ws = walk_size;

   // frame for walk_size loop
   {
      PC_Frame_Base loop("loop", ws);

      {// frame for random int_in_range
        PC_Frame<random_int_cost_f> frame((void*) prng_get_random_int_in_range);
        frame.costf()((NULL));
    }
      
      // frame for Ball_intersectCoord
      {
         PC_Frame<Ball_intersectCoord_cached_cost_f> frame((void*) Ball_intersectCoord_cached);
         frame.costf()(n);
      }
      
      // body intersectCoord
      for (int c = 0; c < bcount; c++) {
         PC_Frame<intersectCoord_cost_f> frame((void*) type[c]->intersectCoord);
         frame.costf()(body[c]);
      }
      pc_stack().log(2*bcount,0,"Update min/max intersection point for last body.");

      {// frame for random double_in_range
         PC_Frame<random_double_in_range_cost_f> frame((void*) prng_get_random_double_in_range);
         frame.costf()(NULL);
      }

      // Reading and writing x[dd] with one add in between
      pc_stack().log(1, 2*sizeof(FT), "x[dd] += t;");
      
      pc_stack().log(5, 2*sizeof(FT), "update squaredNorm for x");
      
      // body intersectCoord
      for(int c = 0; c < bcount; c++) {
         PC_Frame<cacheUpdateCoord_cost_f> frame((void*) type[c]->cacheUpdateCoord);
         frame.costf()(body[c]);
      }
   }
   
}

void walkCoord_coord_4_cost_ref(const int n, int bcount, const void** body, const Body_T** type) {
   int ws = walk_size;

   // frame for walk_size loop
   {
      PC_Frame_Base loop("loop", ws);

      {// frame for random int_in_range
        PC_Frame<random_int_cost_f> frame((void*) prng_get_random_int_in_range);
        frame.costf()((NULL));
    }
      
      // frame for Ball_intersectCoord
      {
         PC_Frame<Ball_intersectCoord_cached_cost_f> frame((void*) Ball_intersectCoord_cached4);
         frame.costf()(n);
      }
      
      // body intersectCoord
      for (int c = 0; c < bcount; c++) {
         PC_Frame<intersectCoord_cost_f> frame((void*) type[c]->intersectCoord4);
         frame.costf()(body[c]);
      }
      pc_stack().log(4*2*bcount,0,"Update min/max intersection point for last body.");

      {// frame for random double_in_range
         PC_Frame<rand256d_cost_f_t> frame((void*) rand256d_f);
         frame.costf()();
      }
      pc_stack().log(3*4, 0, "range from 0..1");

      // Reading and writing x[dd] with one add in between
      pc_stack().log(4, 4*2*sizeof(FT), "x[dd] += t;");
      
      pc_stack().log(5*4, 4*2*sizeof(FT), "update squaredNorm for x");
      
      // body intersectCoord
      for(int c = 0; c < bcount; c++) {
         PC_Frame<cacheUpdateCoord_cost_f> frame((void*) type[c]->cacheUpdateCoord4);
         frame.costf()(body[c]);
      }
   }
}
void walkCoord_coord_8_cost_ref(const int n, int bcount, const void** body, const Body_T** type) {
   int ws = walk_size;

   // frame for walk_size loop
   {
      PC_Frame_Base loop("loop", ws);

      {// frame for random int_in_range
        PC_Frame<random_int_cost_f> frame((void*) prng_get_random_int_in_range);
        frame.costf()((NULL));
    }
      
      // frame for Ball_intersectCoord
      {
         PC_Frame<Ball_intersectCoord_cached_cost_f> frame((void*) Ball_intersectCoord_cached8);
         frame.costf()(n);
      }
      
      // body intersectCoord
      for (int c = 0; c < bcount; c++) {
         PC_Frame<intersectCoord_cost_f> frame((void*) type[c]->intersectCoord8);
         frame.costf()(body[c]);
      }
      pc_stack().log(8*2*bcount,0,"Update min/max intersection point for last body.");

      {// frame for random double_in_range
         PC_Frame<rand256d_cost_f_t> frame((void*) rand256d_f,2);
         frame.costf()();
      }
      pc_stack().log(3*8, 0, "range from 0..1");

      // Reading and writing x[dd] with one add in between
      pc_stack().log(8, 8*2*sizeof(FT), "x[dd] += t;");
      
      pc_stack().log(5*8, 8*2*sizeof(FT), "update squaredNorm for x");
      
      // body intersectCoord
      for(int c = 0; c < bcount; c++) {
         PC_Frame<cacheUpdateCoord_cost_f> frame((void*) type[c]->cacheUpdateCoord8);
         frame.costf()(body[c]);
      }
   }
}

