#include "global.h"
#include "ellipsoid.h"

#ifndef POLYTOPE_H
#define POLYTOPE_H

typedef struct Polytope {
   FT* A; // size line * m
   // layout (row-wise):
   // a0, a1, ... an, [buffering]
   // a0, a1, ... an, [buffering] 
   // end of line: buffering for 32 byte allignment
   
   FT* b; // size m

   //  A*x <= b
   int line; // size of one row plus buffer for allignment
   int n; // dimensions
   int m; // constraints
} Polytope;



Polytope* Polytope_new(int n, int m);

void Polytope_free(const void* o);

void Polytope_print(const void* o);
bool Polytope_inside_ref(const void* o, const FT* v);

void Polytope_intersect_ref(const void* o, const FT* x, const FT* d, FT* t0, FT* t1);
void Polytope_intersectCoord_ref(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache);
void Polytope_intersectCoord_cached_ref(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache);


int  Polytope_cacheAlloc_ref(const void* o);
void Polytope_cacheReset_ref(const void* o, const FT* x, void* cache);
void Polytope_cacheUpdateCoord_ref(const void* o, const int d, const FT dx, void* cache);
bool Polytope_shallowCutOracle_ref(const void* o, const Ellipsoid* e, FT* v, FT* c);
void Polytope_transform_ref(const void* o_in, void* o_out, const Matrix* L, FT* a, FT beta);

#endif
