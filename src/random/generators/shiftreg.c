#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include "shiftreg.h"

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
    return xorshift32();
}

inline uint64_t sr_random_uint64() {
    return xorshift64();
}
