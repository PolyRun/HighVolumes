#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>

#ifndef HEADER_VOLUMES_H
#define HEADER_VOLUMES_H

typedef double FT;
#define FT_EPS DBL_EPSILON
#define FT_MAX DBL_MAX
#define FT_MIN DBL_MIN

// --------------------------------------------- Vectors

// simple vector product
//   assume no memory allignment!
FT dotProduct(const FT* u, const FT* v, const int n);

// --------------------------------------------- Polytope

typedef struct Polytope Polytope;

struct Polytope {
   FT* data; // size m * (n + 1)
   // layout (row-wise):
   // a0, a1, ... an, b 
   // a0, a1, ... an, b 
   // ...
   // basically always normal vector plus b adjacent
   // this is hopefully good for cache, right?
   
   //  A*x <= b

   int n; // dimensions
   int m; // constraints
};

// Constructor
Polytope* Polytope_new(int n, int m);

// Destroctor
void Polytope_free(Polytope* p);

// Setters:
inline void Polytope_set_a(Polytope* p, int i, int x, FT a) {
   // for constraint i set coefficient for variable x to a
   p->data[i * (p->n+1) + x] = a;
}
inline void Polytope_set_b(Polytope* p, int i, FT b) {
   // for constraint i set b
   p->data[i * (p->n+1) + p->n] = b;
}

// Getters
inline FT Polytope_get_a(const Polytope* p, int i, int x) {
   return p->data[i * (p->n+1) + x];
}
inline FT Polytope_get_b(const Polytope* p, int i) {
   return p->data[i * (p->n+1) + p->n];
}
// get pointer to ai
FT* Polytope_get_aV(const Polytope* p, int i);


// Geometric functions:
bool Polytope_inside(const Polytope* p, const FT* v);
// v: vector of length p->n


// intersect Polytope p with line defined by point x (inside p)
//                                           and direction d
// returns intersections: x+d*t0, x+d*t1
void Polytope_intersect(const Polytope* p, const FT* x, const FT* d, FT* t0, FT* t1);

// --------------------------------------------- Volume estimation

// input:
//   n: dimensions
//   r0, r1: inner and outer radius
//   body: convex body, guaranteed to contain B(0,r0) and be contained by B(0,r1)

FT volumeEstimateNormalizedBody(const int n, const FT r0, const FT r1, const void* body);

#endif // HEADER_VOLUMES_H
