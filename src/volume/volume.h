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

typedef struct Body_T Body_T;
#include "polytope/polytope.h"
#include "ellipsoid/ellipsoid.h"
#include "preprocess/preprocess.h"
#include "body/body.h"

#ifndef HEADER_VOLUMES_H
#define HEADER_VOLUMES_H

//#define DEBUG_MSG
//#define DEBUG
//#define PRINT_T
//#define PRINT_TMI



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
