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

   // ----------------- Ellipsoid
   { // test Ellipsoid_eval and Ellipsoid_normal
      Ellipsoid* e = Ellipsoid_new(3); // simple sphere
      FT n[3];
      {
         FT x[3] = {0,0,0};
	 FT eval = Ellipsoid_eval(e, x);
	 assert(eval == 0);
         Ellipsoid_normal(e, x, n);
	 assert(n[0]==0 && n[1]==0 && n[2]==0);
      }
      {
         FT x[3] = {0,0,1};
	 FT eval = Ellipsoid_eval(e, x);
	 assert(eval == 1.0);
         Ellipsoid_normal(e, x, n);
	 assert(n[0]==0 && n[1]==0 && n[2]==2);
      }
      {
         FT x[3] = {0,2.0,0};
	 FT eval = Ellipsoid_eval(e, x);
	 assert(eval == 4.0);
         Ellipsoid_normal(e, x, n);
	 assert(n[0]==0 && n[1]==4 && n[2]==0);
      } 
      {
         FT x[3] = {3.0,4.0,0};
	 FT eval = Ellipsoid_eval(e, x);
	 assert(eval == 25.0);
         Ellipsoid_normal(e, x, n);
	 assert(n[0]==6 && n[1]==8 && n[2]==0);
      }
      Ellipsoid_free(e);
   }
   { // test Ellipsoid_project
      Ellipsoid* e = Ellipsoid_new(3); // simple sphere
      e->a[0] = 5.0;
      e->a[1] = -10.0;
      e->a[2] = 0.5;
      
      FT* A0 = Ellipsoid_get_Ai(e,0);
      FT* A1 = Ellipsoid_get_Ai(e,1);
      FT* A2 = Ellipsoid_get_Ai(e,2);
      A0[0] = 2.0;
      A1[1] = 1.0;
      A2[2] = 0.5;

      for(int t=0; t<100; t++) {
         FT x[3] = {t % 13,t % 11 - 12, std::fmod(t,0.7)};
         Ellipsoid_project(e, x);
         FT eval = Ellipsoid_eval(e, x);
         assert(std::abs(eval-1) < 0.00001);
      }
   }
   
   for(int t=0; t<100; t++) {
      const int n = 20;
      Ellipsoid* e1 = Ellipsoid_new(n); // simple sphere
      Ellipsoid* e2 = Ellipsoid_new(n); // simple sphere
      for(int i=0; i<n; i++) {
         //e1->a[i] = prng_get_random_double_in_range(-10,10);
         e2->a[i] = prng_get_random_double_in_range(-10,10);
         //FT* A1i = Ellipsoid_get_Ai(e1,i);
         //A1i[i] = prng_get_random_double_in_range(0.1,20);
         FT* A2i = Ellipsoid_get_Ai(e2,i);
         A2i[i] = prng_get_random_double_in_range(0.1,20);
      }
      
      {
         FT x[n];
         for(int i=0; i<n; i++) {
	    x[i] = e2->a[i] + prng_get_random_double_normal();
	 }

         FT n1[n];
         FT n2[n];
	 Ellipsoid_minimize(e2,e1,x); // test
         Ellipsoid_normal(e1, x, n1);
         Ellipsoid_normal(e2, x, n2);
	 FT n1_2 = dotProduct(n1,n1, n);
	 FT dot = dotProduct(n1,n2, n);
	 
         FT step[n];
         for(int i=0;i<n;i++) { step[i] = n2[i] - n1[i] * dot / n1_2;}
	 FT step2 = dotProduct(step,step,n);
	 assert(std::abs(step2) < 0.0001);
      }
      

      Ellipsoid_free(e1);
      Ellipsoid_free(e2);
   }
   // -------------------------------- end tests

   #ifdef NDEBUG
   std::cout<< "WARNING: DEBUG DISABLED!\n";
   #else
   std::cout<< "TESTS COMPLETE.\n";
   #endif
}

