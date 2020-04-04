#include <iostream>

extern "C" { // must be included C stlye
#include "../../src/poly/volume.h"
#include "../../src/poly/preprocess.h"
}

#include "../../src/poly/volume_helper.hpp"

#include "../../src/util/cli.hpp"
#include "../../src/util/cli_functions.hpp"

int main(int argc, char** argv) {
   CLI cli(argc,argv,"test_volume_basics");
   CLIFunctionsVolume cliFun(cli);
  
   cliFun.preParse();
   if (!cli.parse()) {return -1;}
   cliFun.postParse();
   
   // -------------------------------- start tests
 
   // Generate new polytope box, 4 dim, 2 radius
   Polytope* box = Polytope_new_box(4,2);
   
   // -- prototype:
   auto o = dynamic_cast<CLIF_Option<xyz_f_t>*>(cliFun.getOption("xyz_f"));
   for(auto it : o->fmap) {
      assert(it.second(box,0.1,box->n));
      std::cout << "fname: " << it.first << " fcall: " << it.second(box,0.1,4) << std::endl;
   }
   // -- end prototype.


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

   // Test Polytope_intersect:
   {
      FT x[4] = {0.0,0.0,0.0,0.0};
      FT d[4] = {0.1,-1.0,-0.9,0.2};
      FT t0,t1;
      Polytope_intersect(box, x, d, &t0, &t1);
      assert(t0==-2.0 && t1==2.0);
   }
   {
      FT x[4] = {0.0,0.0,1.0,0.0};
      FT d[4] = {0.0,0.0,-1.0,0.0};
      FT t0,t1;
      Polytope_intersect(box, x, d, &t0, &t1);
      assert(t0==-1.0 && t1==3.0);
   }
   {
      FT x[4] = {0.0,0.0,1.0,-1.5};
      FT d[4] = {0.0,0.1,0.1,0.5};
      FT t0,t1;
      Polytope_intersect(box, x, d, &t0, &t1);
      assert(t0==-1.0 && t1==7.0);
   }
   {
      FT x[4] = {1.5,1.5,0.0,0.0};
      FT d[4] = {1.0,-1.0,0.0,0.0};
      FT t0,t1;
      Polytope_intersect(box, x, d, &t0, &t1);
      assert(t0==-0.5 && t1==0.5);
   }
   
   Polytope_free(box);
   
   // Check ball volume:
   assert(std::abs(Ball_volume(3,1.0) - 4.189) <= 0.01);
   assert(std::abs(Ball_volume(4,1.0) - 4.935) <= 0.01);
   assert(std::abs(Ball_volume(10,1.0) - 2.550) <= 0.01);
   assert(std::abs(Ball_volume(11,1.0) - 1.884) <= 0.01);


   // -------------------------------- end tests

   #ifdef NDEBUG
   std::cout<< "WARNING: DEBUG DISABLED!\n";
   #else
   std::cout<< "TESTS COMPLETE.\n";
   #endif
}

