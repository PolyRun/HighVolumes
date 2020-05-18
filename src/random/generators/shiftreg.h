#include <stdlib.h>

#ifndef HEADER_PRNG_SHIFTREG_H
#define HEADER_PRNG_SHIFTREG_H

#include <stdlib.h>

void sr_init(void *seed);

uint32_t sr_random_uint32();

uint64_t sr_random_uint64();


#endif // HEADER_PRNG_SHIFTREG_H
