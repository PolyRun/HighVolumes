#include <stdlib.h>
#include <time.h>
#include <float.h>
#include <math.h>
#include "prng.h"

rand_init_f_t rand_init_f = std_init;
rand_f_t rand_f = std_rand;

int rand_chunk_size = 1024;

/**
 * \brief Initializes the prng
 **/
void prng_init(){
    rand_init_f(NULL);
    //srand((unsigned) time(NULL));
}

/**
 * \brief Returns a new random double
 **/
double prng_get_random_double(){
    return ((double) rand_f() / (RAND_MAX)) * (DBL_MAX);
}

/**
 * \brief Returns a new random double in range [0,1)
 **/
double prng_get_random_double_0_1(){
    return ((double) rand_f() / (RAND_MAX));
}

/**
 * \brief Returns a new random double from normal distribution
 **/
double prng_get_random_double_normal() {
    // box-muller method:
    const double u = prng_get_random_double_0_1();
    const double v = prng_get_random_double_0_1();
    const double lnu = log(u);
    const double twopiv = 2.0*M_PI*v;
    const double x = sqrt(-2.0*lnu) * cos(twopiv);
    // const double y = sqrt(-2.0*lnu) * sin(twopiv); // wasted
    return x;
}

/**
 * \brief Returns a new random double in  range [lower_bound, upper_bound)
 * \param lower_bound Lower bound on the value of the random double
 * \param upper_bound Upper bound on the value of the random double
 **/
double prng_get_random_double_in_range(double lower_bound, double upper_bound){
    return ((double) rand_f() / (RAND_MAX)) * (upper_bound-lower_bound) + lower_bound;
}

/**
 * \brief Returns a new random integer
 **/
int prng_get_random_int(){
    return rand_f();
}

/**
 * \brief Returns a new random integer in range [lower_bound, upper_bound]
 * \param lower_bound Lower bound on the value of the random integer
 * \param upper_bound Upper bound on the value of the random integer
 **/
int prng_get_random_int_in_range(int lower_bound, int upper_bound){
    return (rand_f() % (upper_bound-lower_bound+1)) + lower_bound;
}
