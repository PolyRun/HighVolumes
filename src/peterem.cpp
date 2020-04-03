#include <iostream>

extern "C" { // must be included C stlye
#include "poly/volume.h"
#include "poly/cube.h"
#include "beta_cut.h"
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
   std::cout << "ok.\n\n";

   std::cout << "Test Polytope_intersect:\n";
   {
      FT x[4] = {0.0,0.0,0.0,0.0};
      FT d[4] = {0.1,-1.0,-0.9,0.2};
      FT t0,t1;
      Polytope_intersect(box, x, d, &t0, &t1);
      std::cout << t0 << " " << t1 << std::endl;
      assert(t0==-2.0 && t1==2.0);
   }
   {
      FT x[4] = {0.0,0.0,1.0,0.0};
      FT d[4] = {0.0,0.0,-1.0,0.0};
      FT t0,t1;
      Polytope_intersect(box, x, d, &t0, &t1);
      std::cout << t0 << " " << t1 << std::endl;
      assert(t0==-1.0 && t1==3.0);
   }
   {
      FT x[4] = {0.0,0.0,1.0,-1.5};
      FT d[4] = {0.0,0.1,0.1,0.5};
      FT t0,t1;
      Polytope_intersect(box, x, d, &t0, &t1);
      std::cout << t0 << " " << t1 << std::endl;
      assert(t0==-1.0 && t1==7.0);
   }
   {
      FT x[4] = {1.5,1.5,0.0,0.0};
      FT d[4] = {1.0,-1.0,0.0,0.0};
      FT t0,t1;
      Polytope_intersect(box, x, d, &t0, &t1);
      std::cout << t0 << " " << t1 << std::endl;
      assert(t0==-0.5 && t1==0.5);
   }
   std::cout << "ok.\n\n";
   
   {
      std::cout << "Print ball volumes, r=1:\n";
      for(int i=1; i<13; i++) {
         std::cout << "  n=" << i << " " << Ball_volume(i,1.0) << "\n";
      }
   }

   {
      std::cout << "\nEstimate volume of 4dim r=2 box. r0=1.5, r1 = 5:\n";
      FT v = volumeEstimateNormalizedBody(4,1.5,5.0,box);
      std::cout << "estimated volume of box: " << v << " (exact: 256).\n";
   }

   Polytope_free(p);
   std::cout << "ok.\n\n";
   

   Polytope *P;
   polycube(10, &P);
   std::cout << *P << std::endl;
   Polytope_free(P);
   

   #ifdef NDEBUG
   std::cout<< "## WARNING: DEBUG DISABLED!\n";
   #else
   std::cout<< "## TESTS COMPLETE.\n";
   #endif
}

