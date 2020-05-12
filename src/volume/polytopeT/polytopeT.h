typedef struct PolytopeT PolytopeT;

#ifndef POLYTOPET_H
#define POLYTOPET_H

#include "../volume.h"
#include "intersectCoord/intersectCoord.h"

extern Body_T PolytopeT_T;

typedef struct PolytopeT {
   FT* A; // size line * m
   // layout (column-wise):
   // a0, a0, a0, ... [buffering]
   // a1, a1, a1, ... [buffering]
   // .., .., .., ... [buffering]
   // end of line: buffering for 32 byte allignment
   
   FT* b; // size m

   //  A*x <= b
   int line; // size of one column plus buffer for allignment
   int n; // dimensions
   int m; // constraints
} PolytopeT;


PolytopeT* PolytopeT_new(int n, int m);

void PolytopeT_free(const void* o);
void* PolytopeT_clone(const void* o);

void PolytopeT_print(const void* o);
bool PolytopeT_inside_ref(const void* o, const FT* v);

void PolytopeT_intersect_ref(const void* o, const FT* x, const FT* d, FT* t0, FT* t1);
void PolytopeT_intersectCoord_ref(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache);
#include "intersectCoord/intersectCoord.h"// see variants

int  PolytopeT_cacheAlloc_ref(const void* o);
void PolytopeT_cacheReset_ref(const void* o, const FT* x, void* cache);
void PolytopeT_cacheUpdateCoord_ref(const void* o, const int d, const FT dx, void* cache);
bool PolytopeT_shallowCutOracle_ref(const void* o, const Ellipsoid* e, FT* v, FT* c);
void PolytopeT_transform_ref(const void* o_in, void* o_out, const Matrix* L, const FT* a, const FT beta);
void PolytopeT_bounding_ref(const void *B, FT *R2, FT *ori);

// --------------- inline Accessors:
static inline void PolytopeT_set_a(PolytopeT* p, int i, int x, FT a) __attribute__((always_inline));
static inline void PolytopeT_set_a(PolytopeT* p, int i, int x, FT a) {
   p->A[i + (p->line) * x] = a;
}
static inline void PolytopeT_set_b(PolytopeT* p, int i, FT b) __attribute__((always_inline));
static inline void PolytopeT_set_b(PolytopeT* p, int i, FT b) {
   p->b[i] = b;
}

static inline FT PolytopeT_get_a(const PolytopeT* p, int i, int x) __attribute__((always_inline));
static inline FT PolytopeT_get_a(const PolytopeT* p, int i, int x) {
   return p->A[i + (p->line) * x];
}
static inline FT PolytopeT_get_b(const PolytopeT* p, int i) __attribute__((always_inline));
static inline FT PolytopeT_get_b(const PolytopeT* p, int i) {
   return p->b[i];
}



#endif
