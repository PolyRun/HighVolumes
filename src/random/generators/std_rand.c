#include <stdlib.h>
#include <time.h>
#include "std_rand.h"
#include <stdio.h>

static uint32_t *chunk;
static int chunk_ptr;

void std_init(void *seed){
    printf("std_rand init\n");
    srand((unsigned) time(seed));
}

inline uint32_t std_rand(){
    return rand();
}

inline void refill_chunk() {
    for(int i = 0; i < rand_chunk_size; ++i) {
        chunk[i] = rand();
    }
}

void std_init_chunked(void *seed) {
    printf("std_rand_chunked init\n");
    srand((unsigned) time(seed));
    chunk = (uint32_t*) malloc(rand_chunk_size*sizeof(uint32_t));
    refill_chunk();
    chunk_ptr = -1;
}

uint32_t std_rand_chunked() {
    if (chunk_ptr >= rand_chunk_size-1){
        refill_chunk();
        chunk_ptr = -1;
    }
    chunk_ptr++;
    return chunk[chunk_ptr];
}

__m256i std_rand256i() {
   return _mm256_set_epi32(
		   rand(),
		   rand(),
		   rand(),
		   rand(),
		   rand(),
		   rand(),
		   rand(),
		   rand()
		   );
}

__m256d std_rand256d() {
   __m256d rr = _mm256_set_pd(rand(), rand(), rand(), rand());
   const __m256d rMaxInv = _mm256_set1_pd(1.0/RAND_MAX);
   __m256d r = _mm256_mul_pd(rr,rMaxInv);
   return r;
}


