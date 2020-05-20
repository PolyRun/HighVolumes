#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <immintrin.h>

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
    sr_refill_chunk();
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

/* Vectorized */

uint32_t *chunk_vec;
int chunk_ptr_vec;

static __m256i state_32_vec;
static __m256i const_05_vec;
static __m256i const_13_vec;
static __m256i const_17_vec;
static __m256i clear_mask_vec;
static __m256i store_mask_vec;

/* Debug function */
void print_vec(__m256i vec) {
    uint32_t out[8];
    _mm256_maskstore_epi32 (out, store_mask_vec, vec);
    for (int i = 0; i < 8; ++i) {
        printf("Elem %i is %i\n", i, out[i]);
    }
}

void sr_init_vec(void *seed) {
    // Ignore seed for simplicity
    srand((unsigned) time(NULL));
    int r[8];
    for(int i = 0; i < 8; ++i) {
        do {
            r[i] = rand();
        } while (r[i] == 0);
    }
    int mask = (1UL << 31);
    int clear_mask = ~mask;
    state_32_vec = _mm256_set_epi32(r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7]);
    const_05_vec = _mm256_set1_epi32(5);
    const_13_vec = _mm256_set1_epi32(13);
    const_17_vec = _mm256_set1_epi32(17);
    clear_mask_vec = _mm256_set1_epi32(clear_mask);
    store_mask_vec = _mm256_set1_epi32(mask);
    chunk_vec = (uint32_t*) malloc(8*sizeof(uint32_t));
    sr_refill_vec();
    chunk_ptr_vec = -1;
}


void sr_refill_vec() {
    __m256i tmp = _mm256_sllv_epi32(state_32_vec, const_13_vec);
    state_32_vec = _mm256_xor_si256(state_32_vec, tmp);
    tmp = _mm256_srav_epi32(state_32_vec, const_17_vec);
    state_32_vec = _mm256_xor_si256(state_32_vec, tmp);
    tmp = _mm256_sllv_epi32(state_32_vec, const_05_vec);
    state_32_vec = _mm256_xor_si256(state_32_vec, tmp);
    tmp = _mm256_and_si256(state_32_vec, clear_mask_vec);
    _mm256_maskstore_epi32 (chunk_vec, store_mask_vec, tmp);
}

uint32_t sr_rand_vec() {
    if (chunk_ptr_vec >= 7){
        sr_refill_vec();
        chunk_ptr_vec = -1;
    }
    chunk_ptr_vec++;
    return chunk_vec[chunk_ptr_vec];
}
