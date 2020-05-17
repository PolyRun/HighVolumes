#include "../ellipsoid/ellipsoid.h"
#include "../polytopeT/polytopeT.h"
#include "../polytopeCSC/polytopeCSC.h"
#include "../polytopeJIT/polytopeJIT.h"
#include "../ft.h"
#include "../linalg/linalg.h"

#ifndef POLYTOPE_H
#define POLYTOPE_H

extern Body_T Polytope_T;

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
void* Polytope_clone(const void* o);
PolytopeT* Polytope_to_PolytopeT(const Polytope* p);
PolytopeCSC *Polytope_to_PolytopeCSC(const Polytope *p);

void Polytope_print(const void* o);
bool Polytope_inside_ref(const void* o, const FT* v);

void Polytope_intersect_ref(const void* o, const FT* x, const FT* d, FT* t0, FT* t1);
void Polytope_intersectCoord_ref(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache);
void Polytope_intersectCoord_cached_ref(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache);


int  Polytope_cacheAlloc_ref(const void* o);
void Polytope_cacheReset_ref(const void* o, const FT* x, void* cache);
void Polytope_cacheUpdateCoord_ref(const void* o, const int d, const FT dx, void* cache);
bool Polytope_shallowCutOracle_ref(const void* o, const Ellipsoid* e, FT* v, FT* c);
void Polytope_transform_ref(const void* o_in, void* o_out, const Matrix* L, const FT* a, const FT beta);

void Polytope_bounding_ref(const void *B, FT *R2, FT *ori);

void Polytope_normal(const void* o, const FT* x, FT* normal);

// --------------- inline Accessors:
static inline FT* Polytope_get_Ai(const Polytope* p, int i) __attribute__((always_inline));
static inline FT* Polytope_get_Ai(const Polytope* p, int i) {
   //return &(p->A[i * (p->line)]);
   return p->A + (i * (p->line));
}

static inline void Polytope_set_a(Polytope* p, int i, int x, FT a) __attribute__((always_inline));
static inline void Polytope_set_a(Polytope* p, int i, int x, FT a) {
   p->A[i * (p->line) + x] = a;
}
static inline void Polytope_set_b(Polytope* p, int i, FT b) __attribute__((always_inline));
static inline void Polytope_set_b(Polytope* p, int i, FT b) {
   p->b[i] = b;
}

static inline FT Polytope_get_a(const Polytope* p, int i, int x) __attribute__((always_inline));
static inline FT Polytope_get_a(const Polytope* p, int i, int x) {
   return p->A[i * (p->line) + x];
}
static inline FT Polytope_get_b(const Polytope* p, int i) __attribute__((always_inline));
static inline FT Polytope_get_b(const Polytope* p, int i) {
   return p->b[i];
}



#endif
