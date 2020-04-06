#include <iostream>

extern "C" { // must be included C stlye
#include "../../src/poly/volume.h"
#include "../../src/poly/preprocess.h"
}

#include "../../src/poly/volume_helper.hpp"

#include "../../src/util/cli.hpp"
#include "../../src/util/cli_functions.hpp"

int main(int argc, char** argv) {
   CLI cli(argc,argv,"test_volume_estimate");
   CLIFunctionsVolume cliFun(cli);
  
   cliFun.preParse();
   if (!cli.parse()) {return -1;}
   cliFun.postParse();
   
   // -------------------------------- start tests
 
   Polytope* box = Polytope_new_box(4,2);
   
   Polytope_T.print(box);

   //// -- prototype:
   //auto o = dynamic_cast<CLIF_Option<xyz_f_t>*>(cliFun.getOption("xyz_f"));
   //for(auto it : o->fmap) {
   //   assert(it.second(box,0.1,box->n) == box->n);
   //   std::cout << "fname: " << it.first << " fcall: " << it.second(box,0.1,4) << std::endl;
   //}
   //// -- end prototype.


   Polytope_free(box);
   

   // -------------------------------- end tests

   #ifdef NDEBUG
   std::cout<< "WARNING: DEBUG DISABLED!\n";
   #else
   std::cout<< "TESTS COMPLETE.\n";
   #endif
}

