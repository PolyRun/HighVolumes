#include <iostream>

extern "C" { // must be included C stlye
#include "poly/volume.h"
    #include "poly/preprocess.h"
}

#include "peterem.hpp"

#include "util/timer.hpp"
#include "util/cli.hpp"
#include "util/cli_functions.hpp"

int main(int argc, char** argv) {
   Timer t;
   t.start();
   t.stop();
   std::cout << "T0: " << t.millisecs() <<"\n";
   
   CLI cli(argc,argv,"peterem");
   cli.addFlag('v',false,"verbose");
   cli.addOption('n',"100","size of data");
   CLIParameters clip;
   clip.set("c","asdfasdf");
   cli.addParameters('p',clip,"somee extra parameters");

   CLIFunctionsVolume cliFun(cli);
  
   cliFun.preParse();
   if (!cli.parse()) {return -1;}
   cliFun.postParse();

   std::cout << "v flag: " << cli.flag('v') << ", n: " << cli.option('n') << std::endl;
   std::cout << "Print parameters p:\n";
   std::cout << cli.parameters('p').tostring() << std::endl;
   std::cout << "And for c: " << cli.parameters('p').get("c","") << "\n\n";
   
   // ---
   {
      Polytope* p = Polytope_new_box(4,2);
      std::cout << "f1: " << xyz_f(p,0.1,4) << std::endl;
      xyz_f = xyz_f2;
      std::cout << "f2: " << xyz_f(p,0.1,4) << std::endl;
      Polytope_free(p);
   }

   // ----- end proof of concept

   std::cout << "\nGenerate new polytope:\n";
   Polytope* p = Polytope_new(10,6);
   hello(p);

   std::cout << "n: " << p->n << " m: " << p->m << " line: " << p->line << "\n";
   
   Polytope_free(p);


   #ifdef NDEBUG
   std::cout<< "## WARNING: DEBUG DISABLED!\n";
   #else
   std::cout<< "## TESTS COMPLETE.\n";
   #endif
}

