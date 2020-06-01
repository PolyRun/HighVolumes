typedef struct PolytopeCSC PolytopeCSC;

#ifndef POLYTOPECSC_H
#define POLYTOPECSC_H

#include "../volume.h"
#include "intersectCoord/intersectCoord.h"

// CSC compressed sparse column format
typedef struct PolytopeCSC {

    // vector holding values
    FT *A;
    FT *Ainv;
    FT *b;

    // A[i] is in row row_idx[i]
    int *row_idx;
    // A[col_start[i]] is first value of col i iff col_start[i+1] != col_start[i]
    // col_start[n] = #nonzeros in A
    int *col_start;
    int n;
    int m;
    
} PolytopeCSC;


extern Body_T PolytopeCSC_T;




extern FT *dotproduct_store_d;
extern FT *dotproduct_store_x;

int nonzerosCSC(const PolytopeCSC *p);


// WARNING: costly, O(nm) worst case cost
FT *PolytopeCSC_get_Ai(const PolytopeCSC *p, int row, FT *res);

void PolytopeCSC_mvm(const PolytopeCSC *p, const FT *x, FT *res);

void PolytopeCSC_print(const void *o);
void PolytopeCSC_free(const void *o);
void *PolytopeCSC_clone(const void *o);
bool PolytopeCSC_inside_ref(const void *o, const FT *x);
void PolytopeCSC_intersect_ref(const void *o, const FT *x, const FT *dir, FT *t0, FT *t1);
void PolytopeCSC_intersectCoord_ref(const void *o, const FT *x, const int d, FT *t0, FT *t1, void *cache);

int PolytopeCSC_cacheAlloc_ref(const void *o);
void PolytopeCSC_cacheReset_ref(const void *o, const FT *x, void *cache);
void PolytopeCSC_cacheReset_withb(const void *o, const FT *x, void *cache);
void PolytopeCSC_cacheReset_fma(const void *o, const FT *x, void *cache);
void PolytopeCSC_cacheReset_vec(const void *o, const FT *x, void *cache);
void PolytopeCSC_cacheUpdateCoord_ref(const void *o, const int d, const FT dx, void *cache);
void PolytopeCSC_cacheUpdateCoord_withb(const void *o, const int d, const FT dx, void *cache);
void PolytopeCSC_cacheUpdateCoord_fma(const void *o, const int d, const FT dx, void *cache);
void PolytopeCSC_cacheUpdateCoord_vec(const void *o, const int d, const FT dx, void *cache);

bool PolytopeCSC_shallowCutOracle_ref(const void* o, const Ellipsoid* e, FT* v, FT* c);

void PolytopeCSC_transform_ref(const void *o_in, void *o_out, const Matrix *L, const FT *a, const FT beta);

void PolytopeCSC_bounding_ref(const void *o, FT *r, FT *ori);


FTpair4 PolytopeCSC_intersectCoord4_ref(const void* o, const FT* x, const int d, void* cache);
FTpair8 PolytopeCSC_intersectCoord8_ref(const void* o, const FT* x, const int d, void* cache);
void PolytopeCSC_cacheReset4_ref(const void* o, const FT* x, void* cache);
void PolytopeCSC_cacheReset8_ref(const void *p, const FT *x, void *cache);
void PolytopeCSC_cacheUpdateCoord4_ref(const void* o, const int d, const __m256d dx, void* cache);
void PolytopeCSC_cacheUpdateCoord8_ref(const void* o, const int d, const FTset8 dx, void* cache);

#endif
