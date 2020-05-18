#include <stdint.h>
#include <stdlib.h>
#include "shiftreg.h"

uint32_t state_32 = 1;
uint64_t state_64 = 1;

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

void sr_init(int seed){
    state32 = (xorshift32_state*) malloc(sizeof(xorshift32_state));
    state64 = (xorshift64_state*) malloc(sizeof(xorshift64_state));
    if (seed != NULL) {
        state_32 = seed;
        state_64 = seed;
    } else { 
        srand((unsigned) time(seed));
        state_32 = rand();
        state_64 = rand();
    }
}

uint32_t sr_random_uint32(){
    return xorshift32();
}

uint64_t sr_random_uint64(){
    return xorshift64();
}
