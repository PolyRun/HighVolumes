#include <stdlib.h>
#include <stdint.h>

#ifndef HEADER_PRNG_SHIFTREG_H
#define HEADER_PRNG_SHIFTREG_H

extern int rand_chunk_size;

void sr_init(void *seed);

uint32_t sr_random_uint32();

uint64_t sr_random_uint64();


#endif // HEADER_PRNG_SHIFTREG_H
