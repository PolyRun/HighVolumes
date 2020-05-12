#ifndef POLYTOPE_JIT_H
#define POLYTOPE_JIT_H

#include "../volume.h"

extern Body_T PolytopeJIT_T;

// position x
typedef bool (*pjit_inside_f_t)(const FT*);

// direction d, outputs t0,t1
typedef void (*pjit_intersect_f_t)(const FT* x, const FT* d, FT* t0, FT* t1);

// input: x, cache
typedef void (*pjit_cacheReset_f_t)(const FT* x, const FT* cache);

// direction d, outputs t0,t1
typedef void (*pjit_intersectCoord_f_t)(const int d, FT* t0, FT* t1, void* cache);

typedef struct PolytopeJIT {
   pjit_inside_f_t inside;
   pjit_intersect_f_t intersect;
   pjit_cacheReset_f_t cacheReset;
   pjit_intersectCoord_f_t intersectCoord;
   
   int n; // dimensions
   int m; // constraints
} PolytopeJIT;



PolytopeJIT* PolytopeJIT_new(int n, int m);
PolytopeJIT *Polytope_to_PolytopeJIT(const Polytope *p);
#include "jit/inside.h"
#include "jit/intersect.h"
#include "jit/cacheReset.h"
#include "jit/intersectCoord.h"

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
