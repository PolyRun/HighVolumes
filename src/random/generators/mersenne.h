#include <stdlib.h>

#ifndef HEADER_PRNG_MERSENNE_H
#define HEADER_PRNG_MERSENNE_H

uint32_t mt_rand();

void mt_init(void *seed);

#endif // HEADER_PRNG_MERSENNE_H
