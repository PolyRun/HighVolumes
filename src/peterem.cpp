#include <iostream>

extern "C" { // must be included C stlye
#include "poly/volume.h"
#include "poly/cube.h"
#include "beta_cut.h"
}

#include "peterem.hpp"

int main() {
   std::cout << "Generate new polytope:\n";
   Polytope* p = Polytope_new(3,6);
   hello(p);

   std::cout << "Generate new polytope box, 4 dim, 2 radius:\n";
   Polytope* box = Polytope_new_box(4,2);
   hello(box);
   
   std::cout << "Test Polytope_inside:\n";
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

   {
      FT x[4] = {0.0,0.0,0.0,0.0};
      FT d[4] = {0.1,-1.0,-1.1,0.2};
      FT t0,t1;

      Polytope_intersect(box, x, d, &t0, &t1);

      std::cout << "intersect box and line t0: " << t0 << ", t1=" << t1 << "\n";
   }
   
   {
      std::cout << "Test ball volume:\n";
      for(int i=1; i<13; i++) {
         std::cout << "  n=" << i << " " << Ball_volume(i,1.0) << "\n";
      }
   }

   {
      FT v = volumeEstimateNormalizedBody(4,1.0,5.0,box);
      std::cout << "estimated volume of box: " << v << " (exact: 256).\n";
   }

   Polytope_free(p);
   

   std::cout << "\n-------------- Test initEllipsoid:\n";

   int n = 10;
   
   Polytope *P;
   cube(n, &P);

   std::cout << P << "\n";
   
   FT r2;
   FT *ori;
   
   initEllipsoid(P,&r2,&ori);

   printf("R2 is %f\nOri is ", r2);

   for (int i = 0; i < n; i++){
     printf("%f ", ori[i]);
   }
   printf("\n");

   Polytope_free(P);
   

   #ifdef NDEBUG
   std::cout<< "## WARNING: DEBUG DISABLED!\n";
   #else
   std::cout<< "## TESTS COMPLETE.\n";
   #endif
}

