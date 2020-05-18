#include <stdlib.h>

#ifndef HEADER_PRNG_MERSENNE_H
#define HEADER_PRNG_MERSENNE_H

uint32_t mersenne_twister();

void mersenne_twister_init(void *seed);

#endif // HEADER_PRNG_MERSENNE_H
