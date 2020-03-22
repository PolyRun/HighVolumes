#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#ifndef HEADER_VOLUMES_H
#define HEADER_VOLUMES_H

typedef double FT;

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

// Geometric functions:
bool Polytope_inside(const Polytope* p, const FT* v);
// v: vector of length p->n


#endif // HEADER_VOLUMES_H
