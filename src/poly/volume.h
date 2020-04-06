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
#define DEBUG
//#define PRINT_T
//#define PRINT_TMI
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

// input: body, point x, cooordinate i  -  output t0, t1
typedef void (*intersectCoord_f_t)(const void*,const FT*,const int,FT*,FT*);

typedef struct Body_T Body_T;
struct Body_T {
   print_f_t print;
   free_f_t free;
   inside_f_t inside;
   intersect_f_t intersect;
   intersectCoord_f_t intersectCoord;
};

extern Body_T Polytope_T;
extern Body_T Sphere_T;

// --------------------------------------------- Polytope

typedef struct Polytope Polytope;

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

// Setters:
void Polytope_set_a(Polytope* p, int i, int x, FT a);
// for constraint i set coefficient for variable x to a
void Polytope_set_b(Polytope* p, int i, FT b);
// for constraint i set b

// Getters
FT Polytope_get_a(const Polytope* p, int i, int x);
FT Polytope_get_b(const Polytope* p, int i);

// get pointer to ai
FT* Polytope_get_aV(const Polytope* p, int i);

// --------------------------------------------- Sphere

typedef struct Sphere Sphere;

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



// --------------------------------------------- Volume estimation

// input:
//   n: dimensions
//   r0, r1: inner and outer radius
//   body: convex body, guaranteed to contain B(0,r0) and be contained by B(0,r1)
FT volumeEstimateNormalizedBody(const int n, const FT r0, const FT r1, const Polytope* body);
// Note: for now the body is just a polytope, we could make this more generic later!


// ---------- proof of concept:

typedef FT (*xyz_f_t)(const Polytope*,const FT,const int);
FT xyz_f1(const Polytope* body, const FT r, const int n);
FT xyz_f2(const Polytope* body, const FT r, const int n);

extern xyz_f_t xyz_f;

// --- end proof of concept.



#endif // HEADER_VOLUMES_H
