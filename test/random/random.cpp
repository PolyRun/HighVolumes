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

   auto o = dynamic_cast<CLIF_TrippleOption<rand_init_f_t,rand_f_t,rand256d_f_t>*>(cliFun.getOption("rand_f"));
   for(auto it : o->fmap) {
      rand_init_f     = it.second.first.first;
      rand_f          = it.second.first.second.first;
      rand256d_f      = it.second.first.second.second;

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
   

   // Test 4-way random gen:
   for(auto it : o->fmap) {
      rand_init_f     = it.second.first.first;
      rand_f          = it.second.first.second.first;
      rand256d_f      = it.second.first.second.second;

      //if (!(it.second.second == "standard rand" || it.second.second == "shift register rand")) {
      //   continue;
      //}
      rand_init_f(NULL);

      for(int t=0;t<10;t++) {
         __m256d rd = rand256d_f();
         printf("%lf %lf %lf %lf\n",rd[0],rd[1],rd[2],rd[3]);
      }

      //  const __m256i exp = _mm256_set1_epi64x(1023L << 52);
      //  const __m256i mask = _mm256_set1_epi64x(0xFFFFFFFF);
      //  
      //  __m256i r = rand256i_f();
      //  printf("%lx %lx %lx %lx\n",r[0],r[1],r[2],r[3]);
      //  r = _mm256_and_si256(r,mask);
      //  printf("%lx %lx %lx %lx\n",r[0],r[1],r[2],r[3]);
      //  r = _mm256_slli_epi64(r,21); // 1 lat, 1 tp
      //  printf("%lx %lx %lx %lx\n",r[0],r[1],r[2],r[3]);
      //  r = _mm256_or_si256(r,exp); // 1 lat, 2 or 3 throughput
      //  printf("%lx %lx %lx %lx\n",r[0],r[1],r[2],r[3]);
      //  __m256d rd = _mm256_castsi256_pd(r);
      //  printf("%lf %lf %lf %lf\n",rd[0],rd[1],rd[2],rd[3]);
      //  // myHack.l = (myHack.l << 21) | (1023L << 52);
      //  // _mm256_castsi256_pd
   }

   #ifdef NDEBUG
   std::cout<< "WARNING: DEBUG DISABLED!\n";
   #else
   std::cout<< "TESTS COMPLETE.\n";
   #endif
}
