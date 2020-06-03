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
typedef void (*pjit_intersectCoordX_f_t)(const int d, __m256d* t0, __m256d* t1, void* cache);

// direction d, distance travelled dx, cache
typedef void (*pjit_cacheUpdateCoord_f_t)(const int d, const FT dx, void* cache);
typedef void (*pjit_cacheUpdateCoord4_f_t)(const int d, const __m256d dx, void* cache);
typedef void (*pjit_cacheUpdateCoord8_f_t)(const int d, const __m256d dx0, __m256d dx1, void* cache);

typedef struct PolytopeJIT {
   pjit_inside_f_t inside;
   pjit_intersect_f_t intersect;
   pjit_cacheReset_f_t cacheReset;
   pjit_intersectCoord_f_t intersectCoord;
   pjit_cacheUpdateCoord_f_t cacheUpdateCoord;

   pjit_cacheReset_f_t cacheReset4;
   pjit_cacheReset_f_t cacheReset8;
   pjit_intersectCoordX_f_t intersectCoord4;
   pjit_intersectCoordX_f_t intersectCoord8;
   pjit_cacheUpdateCoord4_f_t cacheUpdateCoord4;
   pjit_cacheUpdateCoord8_f_t cacheUpdateCoord8;

   int n; // dimensions
   int m; // constraints

   // pc-count data:
   int nzA;
   size_t intersectCoord_bytes;
   size_t cacheUpdateCoord_bytes;
} PolytopeJIT;

typedef enum PolytopeJIT_Generator {
   pjit_single_rax,     // load single aij at a time, via rax
   pjit_single_data,    // load single aij at a time, via data table
   pjit_single_data_acc,// possibly more than one acc per min/max
   pjit_double_data,    // load double/two aij at a time, via data table
   pjit_quad_data,      // load quad/four aij at a time, via data table
   pjit_quad_data_acc,  // more than one acc per min/max
} PolytopeJIT_Generator;

extern PolytopeJIT_Generator PolytopeJIT_generator;

PolytopeJIT* PolytopeJIT_new(int n, int m);
PolytopeJIT *Polytope_to_PolytopeJIT(const Polytope *p);
#include "jit/inside.h"
#include "jit/intersect.h"
#include "jit/cacheReset.h"
#include "jit/intersectCoord.h"
#include "jit/cacheUpdateCoord.h"

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


FTpair4 PolytopeJIT_intersectCoord4_ref(const void* o, const FT* x, const int d, void* cache);
FTpair8 PolytopeJIT_intersectCoord8_ref(const void* o, const FT* x, const int d, void* cache);
void PolytopeJIT_cacheReset4_ref(const void* o, const FT* x, void* cache);
void PolytopeJIT_cacheReset8_ref(const void *p, const FT *x, void *cache);
void PolytopeJIT_cacheUpdateCoord4_ref(const void* o, const int d, const __m256d dx, void* cache);
void PolytopeJIT_cacheUpdateCoord8_ref(const void* o, const int d, const FTset8 dx, void* cache);

#endif
