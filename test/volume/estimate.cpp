#include <iostream>

extern "C" { // must be included C stlye
#include "../../src/volume/volume.h"
#include "../../src/volume/preprocess.h"
}

#include "../../src/volume/volume_helper.hpp"

#include "../../src/util/cli.hpp"
#include "../../src/util/cli_functions.hpp"

void test() {
   const int n = 4;
   
   Polytope* box = Polytope_new_box(n,2);
   Polytope* box2 = Polytope_new_box(n,1);
   Polytope_T.print(box);
   
   Ellipsoid* e = Ellipsoid_new(n);
   for(int i=0;i<n;i++) {
      FT* Ai = Ellipsoid_get_Ai(e,i);
      Ai[i] = 1.0/9.0;
   }
   Ellipsoid_T.print(e);

   // set up arrays to put in the sub-bodies:
   void* body[2] = {box, e};
   Body_T* type[2] = {&Polytope_T, &Ellipsoid_T};
   
   void* body2[2] = {box, box2};
   Body_T* type2[2] = {&Polytope_T, &Polytope_T};

   for(int i=0; i<5; i++) {
      FT rx = 0.1*i;
      FT r0 = 1.5 - rx;
      FT r1 = 5.0 + rx;
      FT vbox = volume_ref(n,r0,r1,  1,(const void**)&box,  (const Body_T**)&type);
      FT vs   = volume_ref(n,r0,r1,  1,(const void**)&e,    (const Body_T**)&type[1]);
      FT v    = volume_ref(n,r0,r1,  2,(const void**)&body, (const Body_T**)&type);
      FT vbox2 = volume_ref(n,r0*0.5,r1,  2,(const void**)&body2, (const Body_T**)&type2);
      
      std::cout << "\nbox: " << vbox << " vs 256\n";
      FT vs_ref = Ball_volume(n,3.0);
      std::cout << "sphere: " << vs << " vs " << vs_ref << "\n";
      std::cout << "intersection: " << v << "\n";
      std::cout << "\nbox2: " << vbox2 << " vs 16\n";

      assert(std::abs(vbox - std::pow(4,n)) <= 10); // Not great, want to get more accuracy
      assert(std::abs(vs_ref - vs) <= 10);
      assert(vbox > v && vs > v);
      assert(std::abs(vbox2 - 16) <= 1.0);
   }


   Polytope_T.free(box);
   Ellipsoid_T.free(e);
}

int main(int argc, char** argv) {
   CLI cli(argc,argv,"test_volume_estimate");
   CLIFunctionsVolume cliFun(cli);
  
   cliFun.preParse();
   if (!cli.parse()) {return -1;}
   cliFun.postParse();
   
   // -------------------------------- start tests

   auto o = dynamic_cast<CLIF_Option<walk_f_t>*>(cliFun.getOption("walk_f"));
   for(auto it : o->fmap) {
      std::cout << "## run walk_f: " << it.first << std::endl;
      walk_f = it.second; // test for all walk functions
      test();
      std::cout << "## done." << std::endl;
   }
 
   // -------------------------------- end tests

   #ifdef NDEBUG
   std::cout<< "WARNING: DEBUG DISABLED!\n";
   #else
   std::cout<< "TESTS COMPLETE.\n";
   #endif
}

