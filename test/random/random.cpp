#include <iostream>
#include <cassert>
#include <string>
#include <math.h>

#include "../../src/volume/volume_helper.hpp"

extern "C" { // must be included C stlye
#include "../../src/random/prng.h"
}

#define MAY_VAL (int) pow(2,15)
#define NR_NUMBERS (int) pow(2,20)
#define NR_BUCKETS 64
#define EXPECTED_NR (int) NR_NUMBERS/NR_BUCKETS
#define EPS pow(2, -5)


int main(int argc, char *argv[]) {
   CLI cli(argc,argv,"test");
   CLIFunctionsVolume cliFun(cli);

   cliFun.preParse();
   if (!cli.parse()) {return -1;}
   cliFun.postParse();

   std::cout << "start test.\n";

   auto o = dynamic_cast<CLIF_DoubleOption<rand_init_f_t,rand_f_t>*>(cliFun.getOption("rand_f"));
   for(auto it : o->fmap) {
      rand_init_f = it.second.first.first;
      rand_f      = it.second.first.second;

      if (!(it.second.second == "standard rand" || it.second.second == "shift register rand")) {
         continue;
      }
      std::cout << "Check uniformity of " << it.first << " - " << it.second.second << std::endl;
      
      long *buckets = (long*) malloc(NR_BUCKETS*sizeof(long));
      memset(buckets, 0, NR_BUCKETS*sizeof(long));

      rand_f();

      for(int i = 0; i < NR_NUMBERS; ++i) {
         int number = prng_get_random_int();

         buckets[number % NR_BUCKETS]++;
      }  

      for(int i = 0; i < NR_BUCKETS; ++i) {
         assert((1.0*abs(buckets[i]-EXPECTED_NR)/EXPECTED_NR < EPS) && "Not enough uniformity");
      }

      /*for(int i = 0; i < NR_BUCKETS; ++i) {
         if (i < 10) {
            std::cout << "Bucket 0"<< i << ": Numbers inside: " << buckets[i] << ", numbers expected: " << EXPECTED_NR << ", off by: " << 1.0*abs(buckets[i]-EXPECTED_NR)/EXPECTED_NR << std::endl; 
         } else {
            std::cout << "Bucket "<< i << ": Numbers inside: " << buckets[i] << ", numbers expected: " << EXPECTED_NR << ", off by: " << 1.0*abs(buckets[i]-EXPECTED_NR)/EXPECTED_NR << std::endl; 
         }
      }*/
      free(buckets);
   }


   #ifdef NDEBUG
   std::cout<< "WARNING: DEBUG DISABLED!\n";
   #else
   std::cout<< "TESTS COMPLETE.\n";
   #endif
}
