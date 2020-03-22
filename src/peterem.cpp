#include <iostream>

extern "C" { // must be included C stlye
#include "poly/volume.h"
#include "poly/cube.h"
#include "beta_cut.h"
}

#include "peterem.hpp"

int main() {
  int n = 10;
  
  Polytope *P;
  cube(n, &P);

  Polytope_print(P);
  
  FT r2;
  FT *ori;
  
  initEllipsoid(P,&r2,&ori);

  printf("R2 is %f\nOri is ", r2);

  for (int i = 0; i < n; i++){
    printf("%f ", ori[i]);
  }
  printf("\n");

  Polytope_free(P);
}

