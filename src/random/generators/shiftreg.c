#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include "shiftreg.h"

uint32_t *chunk;
int chunk_ptr;

static uint32_t state_32 = 1;
static uint64_t state_64 = 1;

/* Base implementation taken from: https://en.wikipedia.org/wiki/Xorshift#Example_implementation */
/* The state word must be initialized to non-zero */
/* Algorithm "xor" from p. 4 of Marsaglia, "Xorshift RNGs" */

uint32_t xorshift32() {
	state_32 ^= state_32 << 13;
	state_32 ^= state_32 >> 17;
	state_32 ^= state_32 << 5;
	return state_32;
}

uint64_t xorshift64() {
	state_64 ^= state_64 << 13;
	state_64 ^= state_64 >> 7;
	state_64 ^= state_64 << 17;
	return state_64;
}

void sr_init(void *seed_){
    if (seed_ != NULL) {
        uint32_t seed_32 = *((uint32_t*) seed_);
        uint64_t seed_64 = *((uint64_t*) seed_);
        state_32 = seed_32;
        state_64 = seed_64;
    } else { 
        srand((unsigned) time(seed_));
        state_32 = rand();
        state_64 = rand();
    }
}

inline uint32_t sr_random_uint32() {
    return (xorshift32() & ~(1UL << 31)); // Clearing sign bit
}

inline uint64_t sr_random_uint64() {
    return xorshift64();
}

inline void sr_refill_chunk() {
    for(int i = 0; i < rand_chunk_size; ++i) {
        chunk[i] = sr_random_uint32();
    }
}

void sr_init_chunked(void *seed) {
    if (seed != NULL) {
        uint32_t seed_32 = *((uint32_t*) seed);
        state_32 = seed_32;
    } else { 
        srand((unsigned) time(seed));
        state_32 = rand();
    }
    chunk = (uint32_t*) malloc(rand_chunk_size*sizeof(uint32_t));
    chunk_ptr = -1;
}

uint32_t sr_rand_chunked() {
    if (chunk_ptr >= rand_chunk_size-1){
        sr_refill_chunk();
        chunk_ptr = -1;
    }
    chunk_ptr++;
    return chunk[chunk_ptr];
}
