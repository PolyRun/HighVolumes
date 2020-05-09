#include "../volume.h"
//#include "../ellipsoid/ellipsoid.h"
//#include "../polytopeT/polytopeT.h"
//#include "../polytopeCSC/polytopeCSC.h"
#include "../ft.h"
//#include "../linalg/linalg.h"

#ifndef POLYTOPE_JIT_H
#define POLYTOPE_JIT_H

extern Body_T PolytopeJIT_T;

typedef struct PolytopeJIT {
   // TODO: add function pointers here.

   int n; // dimensions
   int m; // constraints
} PolytopeJIT;



PolytopeJIT* PolytopeJIT_new(int n, int m);

void PolytopeJIT_free(const void* o);
void* PolytopeJIT_clone(const void* o);

void PolytopeJIT_print(const void* o);
bool PolytopeJIT_inside_ref(const void* o, const FT* v);

void PolytopeJIT_intersect_ref(const void* o, const FT* x, const FT* d, FT* t0, FT* t1);
void PolytopeJIT_intersectCoord_ref(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache);


int  PolytopeJIT_cacheAlloc_ref(const void* o);
void PolytopeJIT_cacheReset_ref(const void* o, const FT* x, void* cache);
void PolytopeJIT_cacheUpdateCoord_ref(const void* o, const int d, const FT dx, void* cache);
bool PolytopeJIT_shallowCutOracle_ref(const void* o, const Ellipsoid* e, FT* v, FT* c);
void PolytopeJIT_transform_ref(const void* o_in, void* o_out, const Matrix* L, const FT* a, const FT beta);

void PolytopeJIT_bounding_ref(const void *B, FT *R2, FT *ori);


#endif
