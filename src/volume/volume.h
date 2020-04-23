#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <immintrin.h>

#include <assert.h>

#include "../random/prng.h"
#include "cholesky.h"

#ifndef HEADER_VOLUMES_H
#define HEADER_VOLUMES_H

typedef double FT;
#define FT_EPS DBL_EPSILON
#define FT_MAX DBL_MAX
#define FT_MIN DBL_MIN

//#define DEBUG_MSG
//#define DEBUG
//#define PRINT_T
//#define PRINT_TMI

// --------------------------------------------- Forward Declarations
typedef struct Polytope Polytope;
typedef struct Ellipsoid Ellipsoid;

// --------------------------------------------- Vectors / General

// dotProduct:
//   assume no memory allignment!
typedef FT (*dotProduct_f_t)(const FT* u, const FT* v, const int n);
extern dotProduct_f_t dotProduct;
#include "dotProduct/dotProduct.h"

// intersect line x + t*d
//           with ball(0,r)
// assume x is in ball
// return t0,t1 for intersections
void Ball_intersect(const int n, const FT r, const FT* x, const FT* d, FT* t0, FT* t1);

// calculate volume exactly for n-dim ball with radius r
FT Ball_volume(const int n, const FT r);

// given n elements of b bytes, want to get smallest multiple of 32 bytes that fits this.
int ceil_cache(const int n, const int b);


// --------------------------------------------- Memory aligned matrix

// maybe use this inside Polytope... splittig A and b of polytope might make sense for accessing anyway

typedef struct Matrix {
   FT* data; // size m * n
   // layout (row-wise):
   // a00, a01, ... a0n-1, [buffering]
   // ...
   // am-10, am-11, ... am-1n-1, [buffering] 
   int line; // size of one row + b, plus buffer for allignment
   int n; // dimensions
   int m; // constraints   
    
} Matrix;


Matrix* Matrix_new(int n, int m);
void Matrix_free(const void* o);
FT* Matrix_get_row(const Matrix* m, int i);
void Matrix_set(Matrix* m, int i, int x, FT a);
FT Matrix_get(const Matrix* m, int i, int x);
void Matrix_print(const void* o);

// given L (lower triangle), b, solve for x: Lx = b
// simple forward substitution
void Matrix_L_solve(const Matrix* o, FT* x, const FT* b);

// invert positive-definite symmetric matrix In
// and store in out
void Matrix_invert_pdsym(const Matrix *In, Matrix *Out);

// --------------------------------------------- Sub-body Member functions

// input body
typedef void (*print_f_t)(const void*);

// input body
typedef void (*free_f_t)(const void*);

// input: body, point x
typedef bool (*inside_f_t)(const void*,const FT*);

// input: body, point x, direction d  -  output t0, t1
// intersections: x+d*t0, x+d*t1
typedef void (*intersect_f_t)(const void*,const FT*,const FT*,FT*,FT*);

// input: body, point x, cooordinate i, cache  -  output t0, t1
typedef void (*intersectCoord_f_t)(const void*,const FT*,const int,FT*,FT*,void*);

// input: body - output: number of bites cache required
typedef int (*cacheAlloc_f_t)(const void*);

// input: body, vector x, cache
typedef void (*cacheReset_f_t)(const void*, const FT*, void*);

// input: body, dim d, dx on that dim, cache
typedef void (*cacheUpdateCoord_f_t)(const void*, const int, const FT, void*);

// Separation oracle used for preprocessing:
//   If body inside E( (2n)^-2 * A, a):
//       return false
//   else:
//       return true
//       return a plane (v,c) for cutting much of ellipse E(A,a)
//       such that vT * x <= c for all points in body
//       and (2n)^-2 * vT * A * v <= (c - vT * a)^2
//       (cut/touch inner ellipsoid)
//
// input: body, cost/cage ellipsoid
// output: plane (normal v, const c)
//    x in body: vT * x <= c
typedef bool (*shallowCutOracle_f_t)(const void*, const Ellipsoid*, FT*, FT*);

// Transform body after preprocessing
// intput: body_in, body_out, matrix L, vector a, beta.
typedef void (*transform_f_t)(const void*, void*, const Matrix*, FT*, FT);

// input: body (ellipsoid or polytope)
// output: radius FT *r and center FT **ori
// compute sphere with center ori and radius r that encloses the body
typedef void (*boundingSphere_f_t)(const void *, FT *, FT **);

typedef struct Body_T Body_T;
struct Body_T {
   print_f_t print;
   free_f_t free;
   inside_f_t inside;
   intersect_f_t intersect;
   intersectCoord_f_t intersectCoord;
   cacheAlloc_f_t cacheAlloc;
   cacheReset_f_t cacheReset;
   cacheUpdateCoord_f_t cacheUpdateCoord;
   shallowCutOracle_f_t shallowCutOracle;
   transform_f_t transform;
    boundingSphere_f_t boundingSphere;
};

extern Body_T Polytope_T;
extern Body_T Ellipsoid_T;

// --------------------------------------------- Polytope

struct Polytope {
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
};

// Constructor
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
// --------------------------------------------- Ellipsoid

struct Ellipsoid {
   FT* A; // row-wise
   // (x-a)T * A * (x-a)
   FT* a;
   int line; // size of one row of A, plus buffer for allignment
   int n;
   FT* T; // inverse of A
   // only used for preprocessing
   // only non-NULL if generated with Ellipsoid_new_with_T(n)
};

Ellipsoid* Ellipsoid_new(int n);
Ellipsoid* Ellipsoid_new_with_T(int n);// only use if you need T!

void Ellipsoid_free(const void* o);
void Ellipsoid_print(const void* o);
bool Ellipsoid_inside_ref(const void* o, const FT* v);
void Ellipsoid_intersect_ref(const void* o, const FT* x, const FT* d, FT* t0, FT* t1);
void Ellipsoid_intersectCoord_ref(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache);
int  Ellipsoid_cacheAlloc_ref(const void* o);
void Ellipsoid_cacheReset_ref(const void* o, const FT* x, void* cache);
void Ellipsoid_cacheUpdateCoord_ref(const void* o, const int d, const FT dx, void* cache);
bool Ellipsoid_shallowCutOracle_ref(const void* o, const Ellipsoid* e, FT* v, FT* c);
void Ellipsoid_transform_ref(const void* o_in, void* o_out, const Matrix* L, FT* a, FT beta);

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

FT Ellipsoid_eval(const Ellipsoid* e, const FT* x); // (x-a)T * A * (x-a)
void Ellipsoid_normal(const Ellipsoid* e, const FT* x, FT* n); // 2 * A * (x - a)  ==> not normalized

// internal: push x on surface of e: (x-a)T * A * (x-a) = eFac
void Ellipsoid_project(const Ellipsoid* e, const FT eFac, FT* x);


// Constrained to e, minimize f -> write into x
// e = (A,a), f=(B,b)
// (x-a)T * A * (x-a) = eFac
// min_x (x-b)T * B * (x-b)
//
// takes current x as initialization
void Ellipsoid_minimize(const Ellipsoid* e, const FT eFac, const Ellipsoid* f, FT* x);

// recompute A from T (inverse)
void Ellipsoid_A_from_T(Ellipsoid* e);


// --------------------------------------------- Preprocessing

void preprocess_ref(const int n, const int bcount, const void** body_in, void** body_out, const Body_T** type, FT *det);


// --------------------------------------------- Volume estimation

// call this function at before any other function.
// initializes memory arrays
// set max_n to be at least as large as maximum dimension to ever be used
// set max_b to be at least the number of sub-bodies ever used at one time
void volume_lib_init(const int max_n, const int max_b);

// number of points sampled per ball
extern int step_size;
// number of walk-steps taken for a sample
extern int walk_size;

// walk function: n, rk, bcount, body, type, x, d, cache
typedef void (*walk_f_t)(const int, const FT, int bcount, const void**, const Body_T**, FT*, FT*,void**);
void walk_ref(const int n, const FT rk, int bcount, const void** body, const Body_T** type, FT* x, FT* d, void** cache);
void walkCoord_ref(const int n, const FT rk, int bcount, const void** body, const Body_T** type, FT* x, FT* d, void** cache);

extern walk_f_t walk_f;

// input:
//   n: dimensions
//   r0, r1: inner and outer radius
//   bodies: convex body, guaranteed to contain B(0,r0) and be contained by B(0,r1)
//   last 3 arguments: count of bodies, list of bodies, list of body types (member functions)
typedef FT (*volume_f_t)(const int, const FT, const FT, const int, const void**,const Body_T*);

FT volume_ref(const int n, const FT r0, const FT r1, int bcount, const void** body, const Body_T** type);


// ---------- proof of concept:

typedef FT (*xyz_f_t)(const Polytope*,const FT,const int);
FT xyz_f1(const Polytope* body, const FT r, const int n);
FT xyz_f2(const Polytope* body, const FT r, const int n);

extern xyz_f_t xyz_f;

// --- end proof of concept.



#endif // HEADER_VOLUMES_H
