#include <iostream>

extern "C" { // must be included C stlye
#include "../../src/volume/volume.h"
}

#include "../../src/volume/volume_helper.hpp"

#include "../../src/util/cli.hpp"
#include "../../src/util/cli_functions.hpp"

void test() {
   const int n = 4;
   
   Polytope* box = Polytope_new_box(n,2);
   Polytope* box2 = Polytope_new_box(n,1);
   PolytopeT* box3 = PolytopeT_new_box(n,2);
   Polytope_T.print(box);
   
   Ellipsoid* e = Ellipsoid_new_with_T(n);
   for(int i=0;i<n;i++) {
      FT* Ai = Ellipsoid_get_Ai(e,i);
      Ai[i] = 1.0/9.0;
      Ellipsoid_set_Ta(e,i,i,1.0/9.0);
   }
   Ellipsoid_T.print(e);

   // set up arrays to put in the sub-bodies:
   void* body[2] = {box, e};
   Body_T* type[2] = {&Polytope_T, &Ellipsoid_T};
   
   void* body2[2] = {box, box2};
   Body_T* type2[2] = {&Polytope_T, &Polytope_T};
   
   void* body3[1] = {box3};
   Body_T* type3[1] = {&PolytopeT_T};
   
   step_size = 20000; // hot fix this to be a bit faster
   for(int i=0; i<5; i++) {
      FT rx = 0.1*i;
      FT r0 = 1.5 - rx;
      FT r1 = 5.0 + rx;
      ArbitraryExpNum vbox = volume(n,r0,r1,  1,(const void**)&box,  (const Body_T**)&type);
      ArbitraryExpNum vboxT = volume(n,r0,r1,  1,(const void**)&body3,  (const Body_T**)&type3);
      ArbitraryExpNum vs   = volume(n,r0,r1,  1,(const void**)&e,    (const Body_T**)&type[1]);
      ArbitraryExpNum v    = volume(n,r0,r1,  2,(const void**)&body, (const Body_T**)&type);
      ArbitraryExpNum vbox2 = volume(n,r0*0.5,r1,  2,(const void**)&body2, (const Body_T**)&type2);
      
      std::cout << "\nbox: " << vbox.num << "\n";
      std::cout << "\nboxT: " << vboxT.num << "\n";
      FT vs_ref = Ball_volume(n,3.0);
      std::cout << "sphere: " << vs.num << " vs " << vs_ref << "\n";
      std::cout << "intersection: " << v.num << "\n";
      std::cout << "\nbox2: " << vbox2.num << "\n";

      assert(std::abs(vbox.num - std::pow(4,n)) <= 10); // Not great, want to get more accuracy
      assert(std::abs(vboxT.num - std::pow(4,n)) <= 10); 
      assert(std::abs(vs_ref - vs.num) <= 10);
      assert(vbox.num > v.num && vs.num > v.num);
      assert(std::abs(vbox2.num - std::pow(2,n)) <= 1.0);
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

   auto o = dynamic_cast<CLIF_DoubleOption<walk_f_t,volume_f_t>*>(cliFun.getOption("walk_f"));
   for(auto it : o->fmap) {
      std::cout << "## run walk_f: " << it.first << " - " << it.second.second << std::endl;
      walk_f = it.second.first.first; // test for all walk functions
      volume = it.second.first.second; // test for all walk functions
      if(volume == volume_coord_4) {continue;}
      if(volume == volume_coord_8) {continue;}// these are not implemented everyhwere.
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

