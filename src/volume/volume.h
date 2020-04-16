#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>

#include <assert.h>

#include "../random/prng.h"


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
typedef struct Sphere Sphere;
typedef struct Ellipsoid Ellipsoid;

// --------------------------------------------- Vectors / General

// simple vector product
//   assume no memory allignment!
FT dotProduct(const FT* u, const FT* v, const int n);

// intersect line x + t*d
//           with ball(0,r)
// assume x is in ball
// return t0,t1 for intersections
void Ball_intersect(const int n, const FT r, const FT* x, const FT* d, FT* t0, FT* t1);

// calculate volume exactly for n-dim ball with radius r
FT Ball_volume(const int n, const FT r);

// given n elements of b bytes, want to get smallest multiple of 32 bytes that fits this.
int ceil_cache(const int n, const int b);

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
};

extern Body_T Polytope_T;
extern Body_T Sphere_T;
extern Body_T Ellipsoid_T;

// --------------------------------------------- Polytope

struct Polytope {
   FT* data; // size m * (n + 1)
   // layout (row-wise):
   // a0, a1, ... an, b, [buffering]
   // a0, a1, ... an, b, [buffering] 
   // ...
   // basically always normal vector plus b adjacent
   // this is hopefully good for cache, right?
   // end of line: buffering for 32 byte allignment
   
   //  A*x <= b
   int line; // size of one row + b, plus buffer for allignment
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

// Setters:
void Polytope_set_a(Polytope* p, int i, int x, FT a);
// for constraint i set coefficient for variable x to a
void Polytope_set_b(Polytope* p, int i, FT b);
// for constraint i set b

// Getters
FT Polytope_get_a(const Polytope* p, int i, int x);
FT Polytope_get_b(const Polytope* p, int i);

// get pointer to ai
FT* Polytope_get_Ai(const Polytope* p, int i);

// --------------------------------------------- Sphere

struct Sphere {
   FT* center; // size n
   FT r; // radius
   int n; // dimensions
};

Sphere* Sphere_new(int n, FT r, const FT* c);

void Sphere_free(const void* o);
void Sphere_print(const void* o);
bool Sphere_inside_ref(const void* o, const FT* v);
void Sphere_intersect_ref(const void* o, const FT* x, const FT* d, FT* t0, FT* t1);
void Sphere_intersectCoord_ref(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache);
int  Sphere_cacheAlloc_ref(const void* o);
void Sphere_cacheReset_ref(const void* o, const FT* x, void* cache);
void Sphere_cacheUpdateCoord_ref(const void* o, const int d, const FT dx, void* cache);

// --------------------------------------------- Ellipsoid

struct Ellipsoid {
   FT* A; // row-wise
   // (x-a)T * A * (x-a)
   FT* a;
   int line; // size of one row of A, plus buffer for allignment
   int n;
};

Ellipsoid* Ellipsoid_new(int n);

void Ellipsoid_free(const void* o);
void Ellipsoid_print(const void* o);
bool Ellipsoid_inside_ref(const void* o, const FT* v);
void Ellipsoid_intersect_ref(const void* o, const FT* x, const FT* d, FT* t0, FT* t1);
void Ellipsoid_intersectCoord_ref(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache);
int  Ellipsoid_cacheAlloc_ref(const void* o);
void Ellipsoid_cacheReset_ref(const void* o, const FT* x, void* cache);
void Ellipsoid_cacheUpdateCoord_ref(const void* o, const int d, const FT dx, void* cache);
bool Ellipsoid_shallowCutOracle_ref(const void* o, const Ellipsoid* e, FT* v, FT* c);


FT* Ellipsoid_get_Ai(const Ellipsoid* e, int i); // get row i
FT Ellipsoid_eval(const Ellipsoid* e, const FT* x); // (x-a)T * A * (x-a)
void Ellipsoid_normal(const Ellipsoid* e, const FT* x, FT* n); // 2 * A * (x - a)  ==> not normalized

// internal: push x on surface of e
void Ellipsoid_project(const Ellipsoid* e, FT* x);


// Constrained to e, minimize f -> write into x
// e = (A,a), f=(B,b)
// (x-a)T * A * (x-a) = 1
// min_x (x-b)T * B * (x-b)
//
// takes current x as initialization
void Ellipsoid_minimize(const Ellipsoid* e, const Ellipsoid* f, FT* x);

// --------------------------------------------- Preprocessing

void preprocess_ref(const int n, const int bcount, const void** body_in, void** body_out, const Body_T** type, FT *det);


// --------------------------------------------- Volume estimation

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
