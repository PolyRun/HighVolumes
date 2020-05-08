typedef struct PolytopeCSC PolytopeCSC;

#ifndef POLYTOPECSC_H
#define POLYTOPECSC_H

#include "../volume.h"

// CSC compressed sparse column format
typedef struct PolytopeCSC {

    // vector holding values
    FT *A;
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


// WARNING: costly, O(nm) worst case cost
FT *PolytopeCSC_get_Ai(const PolytopeCSC *p, int row, FT *res);


void PolytopeCSC_print(const void *o);
void PolytopeCSC_free(const void *o);
void *PolytopeCSC_clone(const void *o);
bool PolytopeCSC_inside_ref(const void *o, const FT *x);
void PolytopeCSC_intersect_ref(const void *o, const FT *x, const FT *dir, FT *t0, FT *t1);
void PolytopeCSC_intersectCoord_ref(const void *o, const FT *x, const int d, FT *t0, FT *t1, void *cache);

int PolytopeCSC_cacheAlloc_ref(const void *o);
void PolytopeCSC_cacheReset_ref(const void *o, const FT *x, void *cache);
void PolytopeCSC_cacheUpdateCoord_ref(const void *o, const int d, const FT dx, void *cache);

bool PolytopeCSC_shallowCutOracle_ref(const void* o, const Ellipsoid* e, FT* v, FT* c);

void PolytopeCSC_transform_ref(const void *o_in, void *o_out, const Matrix *L, const FT *a, const FT beta);

void PolytopeCSC_bounding_ref(const void *o, FT *r, FT *ori);


#endif
