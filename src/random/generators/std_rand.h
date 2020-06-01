
#ifndef HEADER_PRNG_STD_H
#define HEADER_PRNG_STD_H
#include <stdint.h>
#include <immintrin.h>

extern int rand_chunk_size;

void std_init(void *seed);

uint32_t std_rand();

// Chunked random values
void std_init_chunked(void *seed);

uint32_t std_rand_chunked();

void refill_chunk();

__m256i std_rand256i();
__m256d std_rand256d();
#endif // HEADER_PRNG_STD_H
