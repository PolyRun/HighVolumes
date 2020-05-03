// Please ignore this file, it is only used for small WIP tests

#include <iostream>

extern "C" { // must be included C stlye
#include "volume/volume.h"
    #include "volume/preprocess/preprocess.h"
}

#include "peterem.hpp"

#include "util/timer.hpp"
#include "util/cli.hpp"
#include "util/cli_functions.hpp"

#include "util/image_bmp.hpp"



int main(int argc, char** argv) {
   Timer t;
   t.start();
   t.stop();
   std::cout << "T0: " << t.millisecs() <<"\n";
   
   CLI cli(argc,argv,"peterem");
   cli.addFlag('v',false,"verbose");
   cli.addOption('s',"100","size of data");
   CLIParameters clip;
   clip.set("c","asdfasdf");
   cli.addParameters('p',clip,"somee extra parameters");

   CLIFunctionsVolume cliFun(cli);
  
   cliFun.preParse();
   if (!cli.parse()) {return -1;}
   cliFun.postParse();

   std::cout << "v flag: " << cli.flag('v') << ", n: " << cli.option('s') << std::endl;
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

   std::cout << "n: " << p->n << " m: " << p->m << " line: " << p->line << "\n";
   
   Polytope_free(p);
   
   std::cout << "ceil_cache: " << ceil_cache(31,1) << "\n";
   std::cout << "ceil_cache: " << ceil_cache(32,1) << "\n";
   std::cout << "ceil_cache: " << ceil_cache(33,1) << "\n";
   std::cout << "ceil_cache: " << ceil_cache(1,sizeof(FT)) << "\n";
   std::cout << "ceil_cache: " << ceil_cache(4,sizeof(FT)) << "\n";
   std::cout << "ceil_cache: " << ceil_cache(5,sizeof(FT)) << "\n";
   
   {// ---------------------- IMG
      int width = 200;
      int height = 200;
      EVP::Image_BMP img(height,width);
      for(int i=0; i<height; i++){
          for(int j=0; j<width; j++){
              int r = (unsigned char)((double)i/height*255); ///red
              int g = (unsigned char)((double)j/width*255); ///green
              int b = (unsigned char)(((double)i+j)/(height+width)*255); ///blue
              img.set(j,i,r,g,b);
	  }
      }
      img.toFile("peterem_out.bmp");
   }// ---------------------- END IMG



   #ifdef NDEBUG
   std::cout<< "## WARNING: DEBUG DISABLED!\n";
   #else
   std::cout<< "## TESTS COMPLETE.\n";
   #endif
}




