#include <stdlib.h>
#include <time.h>
#include <float.h>
#include <math.h>
#include "prng.h"

rand_init_f_t rand_init_f = std_init;
rand_f_t rand_f = std_rand;
rand256i_f_t rand256i_f = std_rand256i;
rand256d_f_t rand256d_f = std_rand256d;

rd_0_1_f_t prng_get_random_double_0_1 = prng_get_random_double_0_1_ref;

int rand_chunk_size = 1024;

union hack {
    long l;
    double d;
};

/**
 * \brief Initializes the prng
 **/
void prng_init(){
    rand_init_f(NULL);
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
double prng_get_random_double_0_1_ref(){
    return ((double) rand_f() / (RAND_MAX));
}

double prng_fast_32_get_random_double_0_1() {
    union hack myHack;
    myHack.l = (long) rand_f();
    double rand_double;

    // We put the 32 bits of num to the beginning (MSB-wise) of the mantissa and
    // the 20 bits in the bottom are zero, which is an error of 2^(-20)
    // Note that we can't convert directly to double and instead do address-magic
    myHack.l = (myHack.l << 21) | (1023L << 52);
    rand_double = myHack.d;

    // Since rand_double has a zero exponent (2^0), it is between 2 and 1.
    rand_double = rand_double - 1;
    return rand_double;
}

/*
double prng_fast_64_get_random_double_0_1() {
    long num = (long) xorshift64();
    long mask = (1L << 52) - 1; // This gives a mask which has 12 bits 0 and 52 bits 1
    double rand_double;

    // Note that we can't convert directly to double and instead do address-magic
    num = num & mask;
    rand_double = *(double*) &num; // This is undefined behaviour btw :)

    // Since rand_double has a zero exponent (2^0), it is between 2 and 1
    rand_double = rand_double - 1;
    return rand_double;
}
*/

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
