#include <iostream>

extern "C" { // must be included C stlye
#include "../../src/volume/volume.h"
#include "../../src/volume/preprocess.h"
}

#include "../../src/volume/volume_helper.hpp"

#include "../../src/util/cli.hpp"
#include "../../src/util/cli_functions.hpp"

int main(int argc, char** argv) {
   CLI cli(argc,argv,"test_volume_basics");
   CLIFunctionsVolume cliFun(cli);
  
   cliFun.preParse();
   if (!cli.parse()) {return -1;}
   cliFun.postParse();
   
   // -------------------------------- start tests
   auto o = dynamic_cast<CLIF_Option<intersectCoord_f_t>*>(cliFun.getOption("Polytope_intersectCoord"));
   for(auto it : o->fmap) {
      Polytope_T.intersectCoord = it.second;
      std::cout << "Test Polytope for intersectCoord " << it.first << std::endl;

      // Generate new polytope box, 4 dim, 2 radius
      Polytope* box = Polytope_new_box(4,2);
      {
         FT v[4] = {0,0,0,0};
         assert(Polytope_T.inside(box, (FT*)&v));
      }
      {
         FT v[4] = {0,0,0,3};
         assert(!Polytope_T.inside(box, (FT*)&v));
      }
      {
         FT v[4] = {2,2,2,2};
         assert(Polytope_T.inside(box, (FT*)&v));
      }
      {
         FT v[4] = {1,-1,1,-1};
         assert(Polytope_T.inside(box, (FT*)&v));
      }
      {
         FT v[4] = {-3,0,0,0};
         assert(Polytope_T.inside(box, (FT*)&v));
      }

      // Test Polytope_T.intersect:
      {
         FT x[4] = {0.0,0.0,0.0,0.0};
         FT d[4] = {0.1,-1.0,-0.9,0.2};
         FT t0,t1;
         Polytope_T.intersect(box, x, d, &t0, &t1);
         assert(t0==-2.0 && t1==2.0);
      }
      {
         FT x[4] = {0.0,0.0,1.0,0.0};
         FT d[4] = {0.0,0.0,-1.0,0.0};
         FT t0,t1;
         Polytope_T.intersect(box, x, d, &t0, &t1);
         assert(t0==-1.0 && t1==3.0);
      }
      {
         FT x[4] = {0.0,0.0,1.0,-1.5};
         FT d[4] = {0.0,0.1,0.1,0.5};
         FT t0,t1;
         Polytope_T.intersect(box, x, d, &t0, &t1);
         assert(t0==-1.0 && t1==7.0);
      }
      {
         FT x[4] = {1.5,1.5,0.0,0.0};
         FT d[4] = {1.0,-1.0,0.0,0.0};
         FT t0,t1;
         Polytope_T.intersect(box, x, d, &t0, &t1);
         assert(t0==-0.5 && t1==0.5);
      }
      // Test Polytope_T.intersectCoord
      void* cache = aligned_alloc(32, Polytope_T.cacheAlloc(box));
      {
         FT x[4] = {0.0,0.0,0.0,0.0};
         Polytope_T.cacheReset(box,x,cache);
         for(int d=0;d<4;d++) {
            FT t0,t1;
            Polytope_T.intersectCoord(box, x, 0, &t0, &t1, cache);
            assert(t0==-2.0 && t1==2.0);
         }
      }
      {
         FT x[4] = {0.0,0.0,1.0,0.0};
         Polytope_T.cacheReset(box,x,cache);
         FT t0,t1;
         Polytope_T.intersectCoord(box, x, 2, &t0, &t1, cache);
         assert(t0==-3.0 && t1==1.0);
         Polytope_T.intersectCoord(box, x, 0, &t0, &t1, cache);
         assert(t0==-2.0 && t1==2.0);
         Polytope_T.intersectCoord(box, x, 1, &t0, &t1, cache);
         assert(t0==-2.0 && t1==2.0);
         Polytope_T.intersectCoord(box, x, 3, &t0, &t1, cache);
         assert(t0==-2.0 && t1==2.0);
      }
      free(cache);
      Polytope_free(box);
   }

   // Check ball volume:
   assert(std::abs(Ball_volume(3,1.0) - 4.189) <= 0.01);
   assert(std::abs(Ball_volume(4,1.0) - 4.935) <= 0.01);
   assert(std::abs(Ball_volume(10,1.0) - 2.550) <= 0.01);
   assert(std::abs(Ball_volume(11,1.0) - 1.884) <= 0.01);

   // -------------- Sphere:
   {   
      FT center[4] = {0,1,0,0};
      Sphere* s = Sphere_new(4,3.0,center);
      //Sphere_T.print(s);
      {
         FT v[4] = {-1,0,0,0};
         assert(Sphere_T.inside(s, (FT*)&v));
      }
      {
         FT v[4] = {-3,0,0,0};
         assert(!Sphere_T.inside(s, (FT*)&v));
      }
      {
         FT v[4] = {0,3.99,0,0};
         assert(Sphere_T.inside(s, (FT*)&v));
      }
      {
         FT v[4] = {0,4.01,0,0};
         assert(!Sphere_T.inside(s, (FT*)&v));
      }

      void* cache = aligned_alloc(32, Sphere_T.cacheAlloc(s));
      {
         FT x[4] = {0.0,1.0,0.0,0.0};
         Sphere_T.cacheReset(s,x,cache);
         FT d[4] = {1.0,0.0,0.0,0.0};
         FT t0,t1;
         Sphere_T.intersect(s, x, d, &t0, &t1);
         assert(t0==-3.0 && t1==3.0);
         for(int dd=0;dd<4;dd++) {
            Sphere_T.intersectCoord(s, x, dd, &t0, &t1, cache);
            assert(t0==-3.0 && t1==3.0);
         }
      }
      {
         FT x[4] = {0.0,0.0,0.0,0.0};
         Sphere_T.cacheReset(s,x,cache);
         FT d[4] = {0.0,1.0,0.0,0.0};
         FT t0,t1;
         Sphere_T.intersect(s, x, d, &t0, &t1);
         assert(t0==-2.0 && t1==4.0);
         Sphere_T.intersectCoord(s, x, 1, &t0, &t1, cache);
         assert(t0==-2.0 && t1==4.0);
      }
      {
         FT x[4] = {0.0,3.0,0.0,0.0};
         Sphere_T.cacheReset(s,x,cache);
         FT d[4] = {0.0,1.0,0.0,0.0};
         FT t0,t1;
         Sphere_T.intersect(s, x, d, &t0, &t1);
         assert(t0==-5.0 && t1==1.0);
         Sphere_T.intersectCoord(s, x, 1, &t0, &t1, cache);
         assert(t0==-5.0 && t1==1.0);
      }
      {
         FT x[4] = {0.0,1.0,0.0,0.0};
         Sphere_T.cacheReset(s,x,cache);
         FT d[4] = {1.0,1.0,1.0,1.0};
         FT t0,t1;
         Sphere_T.intersect(s, x, d, &t0, &t1);
         assert(t0==-1.5 && t1==1.5);
      }
      free(cache);

      Sphere_T.free(s);
   }
   // -------------------------------- end tests

   #ifdef NDEBUG
   std::cout<< "WARNING: DEBUG DISABLED!\n";
   #else
   std::cout<< "TESTS COMPLETE.\n";
   #endif
}

