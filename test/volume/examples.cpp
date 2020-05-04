#include <iostream>

extern "C" { // must be included C stlye
#include "../../src/volume/volume.h"
}

#include "../../src/volume/volume_helper.hpp"

#include "../../src/util/cli.hpp"
#include "../../src/util/cli_functions.hpp"

void test_intersect(Solved_Body* sb) {

   const int n = sb->n;
   FT* x = (FT*)(aligned_alloc(32, n*sizeof(FT)));
   for(int i = 0; i < n; i++) {
      x[i] = 0;
   }

   FT* d = (FT*)(aligned_alloc(32, n*sizeof(FT)));

   for(int t = 0; t < n*200; t++) {
      for(int i = 0; i < n; i++) {
         d[i] = prng_get_random_double_normal();
      }
      FT d2 = squaredNorm(d,n);
      
      FT t0 = -FT_MAX;
      FT t1 = FT_MAX;

      for(int b = 0; b < sb->bcount; b++) {
         FT t0_,t1_;
         sb->type[b]->intersect(sb->body[b], x, d, &t0_, &t1_);
         t0 = std::max(t0,t0_);
         t1 = std::min(t1,t1_);
      }
      FT tmax = std::max(t0*t0, t1*t1)*d2;
      FT tmin = std::min(t0*t0, t1*t1)*d2;
      assert(tmin >= 1.0 - 0.00001);
      assert(tmax <= 4*n*n + 0.00001);
   }
}


int main(int argc, char** argv) {
   CLI cli(argc,argv,"test_volume_basics");
   CLIFunctionsVolume cliFun(cli);
  
   cliFun.preParse();
   if (!cli.parse()) {return -1;}
   cliFun.postParse();
   
   // -------------------------------- start tests
   
   // iterate over all example bodies, generate one
   // if it promises to be normalized, check it
   
   for(auto it : solved_body_generator()->gen()) {
      std::cout << "Check body: " << it.first << "\n";
      auto gen = it.second;
      auto sb = gen();
      if(!sb->is_normalized) {
         std::cout << "not normalized - not tested.\n";
      } else {
	 test_intersect(sb);
      }
   }

   // -------------------------------- end tests

   #ifdef NDEBUG
   std::cout<< "WARNING: DEBUG DISABLED!\n";
   #else
   std::cout<< "TESTS COMPLETE.\n";
   #endif
}

