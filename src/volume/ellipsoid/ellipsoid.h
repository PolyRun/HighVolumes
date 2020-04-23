#include "../ft.h"
#include "../linalg/linalg.h"

#ifndef ELLIPSOID_H
#define ELLIPSOID_H

typedef struct Ellipsoid {
   FT* A; // row-wise
   // (x-a)T * A * (x-a)
   FT* a;
   int line; // size of one row of A, plus buffer for allignment
   int n;
   FT* T; // inverse of A
   // only used for preprocessing
   // only non-NULL if generated with Ellipsoid_new_with_T(n)
} Ellipsoid;


Ellipsoid* Ellipsoid_new(int n);
Ellipsoid* Ellipsoid_new_with_T(int n);
void Ellipsoid_free(const void* o);
void Ellipsoid_print(const void* o);
FT Ellipsoid_eval(const Ellipsoid* e, const FT* x);
void Ellipsoid_normal(const Ellipsoid* e, const FT* x, FT* normal);
void Ellipsoid_project(const Ellipsoid* e, const FT eFac, FT* x);
void Ellipsoid_minimize(const Ellipsoid* e, const FT eFac, const Ellipsoid* f, FT* x);
void Ellipsoid_A_from_T(Ellipsoid* e);


bool Ellipsoid_inside_ref(const void* o, const FT* v);

void Ellipsoid_intersect_ref(const void* o, const FT* x, const FT* d, FT* t0, FT* t1);

void Ellipsoid_intersectCoord_ref(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache);

int  Ellipsoid_cacheAlloc_ref(const void* o);

void Ellipsoid_cacheReset_ref(const void* o, const FT* x, void* cache);

void Ellipsoid_cacheUpdateCoord_ref(const void* o, const int d, const FT dx, void* cache);

bool Ellipsoid_shallowCutOracle_ref(const void* o, const Ellipsoid* e, FT* v, FT* c);

void Ellipsoid_transform_ref(const void* o_in, void* o_out, const Matrix* L, FT* a, FT beta);



// --------------------------------------------- Ellipsoid
// get row i of A
static inline FT* Ellipsoid_get_Ai(const Ellipsoid* e, int i) __attribute__((always_inline));
static inline FT* Ellipsoid_get_Ai(const Ellipsoid* e, int i) {
   return e->A + i*e->line;
}

// get row i of T
static inline FT* Ellipsoid_get_Ti(const Ellipsoid* e, int i) __attribute__((always_inline));
static inline FT* Ellipsoid_get_Ti(const Ellipsoid* e, int i) {
   assert(e->T && "T must be allocated");
   return e->T + i*e->line;
}



#endif
