#include "generators/std_rand.h"
#include "generators/mersenne.h"
#include "generators/shiftreg.h"


#ifndef HEADER_PRNG_H
#define HEADER_PRNG_H

// Generator for randomness
typedef void(*rand_init_f_t)(void *seed);
typedef uint32_t(*rand_f_t)();

extern rand_init_f_t rand_init_f;
extern rand_f_t rand_f;

extern int rand_chunk_size;

/**
 * \brief Initializes the prng
 **/
void prng_init();

/**
 * \brief Returns a new random double
 **/
double prng_get_random_double();

/**
 * \brief Returns a new random double in range [0,1]
 **/
double prng_get_random_double_0_1();

/**
 * \brief Returns a new random double from normal distribution
 **/
double prng_get_random_double_normal();

/**
 * \brief Returns a new random double in  range [lower_bound, upper_bound)
 * \param lower_bound Lower bound on the value of the random double
 * \param upper_bound Upper bound on the value of the random double
 **/
double prng_get_random_double_in_range(double lower_bound, double upper_bound);

/**
 * \brief Returns a new random integer
 **/
int prng_get_random_int();

/**
 * \brief Returns a new random integer in range [lower_bound, upper_bound]
 * \param lower_bound Lower bound on the value of the random integer
 * \param upper_bound Upper bound on the value of the random integer
 **/
int prng_get_random_int_in_range(int lower_bound, int upper_bound);

#endif // HEADER_PRNG_H
