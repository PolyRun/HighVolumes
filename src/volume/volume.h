#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <immintrin.h>

#include <assert.h>

#include "../random/prng.h"
#include "ft.h"
#include "linalg/linalg.h"

#include "matrix/matrix.h"
#include "polytope/polytope.h"
#include "ellipsoid/ellipsoid.h"
#include "preprocess/preprocess.h"

#ifndef HEADER_VOLUMES_H
#define HEADER_VOLUMES_H

//#define DEBUG_MSG
//#define DEBUG
//#define PRINT_T
//#define PRINT_TMI

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
