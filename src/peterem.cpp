<<<<<<< HEAD
#include <iostream>

extern "C" { // must be included C stlye
#include "poly/volume.h"
#include "poly/cube.h"
#include "beta_cut.h"
}

#include "peterem.hpp"

int main() {
   Polytope* p = Polytope_new(3,6);
   hello(p);

   Polytope* box = Polytope_new_box(4,2);
   hello(box);
   
   {
      FT v[4] = {0,0,0,0};
      assert(Polytope_inside(box, (FT*)&v));
   }
   {
      FT v[4] = {0,0,0,3};
      assert(!Polytope_inside(box, (FT*)&v));
   }
   {
      FT v[4] = {2,2,2,2};
      assert(Polytope_inside(box, (FT*)&v));
   }
   {
      FT v[4] = {1,-1,1,-1};
      assert(Polytope_inside(box, (FT*)&v));
   }
   {
      FT v[4] = {-3,0,0,0};
      assert(Polytope_inside(box, (FT*)&v));
   }

   Polytope_free(p);

   assert(false && "Expected assert, all tests passed.");
   std::cout << "WARNING: asserts were disabled, no tests run!" << std::endl;

   
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

