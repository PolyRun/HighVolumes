#include <stdlib.h>
#include <stdint.h>

#ifndef HEADER_PRNG_SHIFTREG_H
#define HEADER_PRNG_SHIFTREG_H

extern int rand_chunk_size;

void sr_init(void *seed);

uint32_t sr_random_uint32();

uint64_t sr_random_uint64();

// Chunked random values
void sr_init_chunked(void *seed);

uint32_t sr_rand_chunked();

void sr_refill_chunk();

// Vectorized random values
void sr_init_vec(void *seed);

uint32_t sr_rand_vec();

void sr_refill_vec();

__m256i sr_rand256i();
__m256d sr_rand256d();

#endif // HEADER_PRNG_SHIFTREG_H
