#include <stdio.h>
#include <stdlib.h>

#ifndef HEADER_VOLUMES_H
#define HEADER_VOLUMES_H

typedef float FT;

typedef struct Polytope Polytope;

struct Polytope {
   FT* data;
   int n; // dimensions
   int m; // constraints
};

Polytope* Polytope_new(int n, int m);

void Polytope_free(Polytope* p);

//Polytope* Polytope_new(int n, int m) {
//   Polytope* o = (Polytope*) malloc(sizeof(Polytope));
//   o->data = (FT*) malloc(sizeof(FT)*(n+1)*m);
//   o->n = n;
//   o->m = m;
//
//   return o;
//}



#endif // HEADER_VOLUMES_H
