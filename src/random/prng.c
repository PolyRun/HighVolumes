#include <stdlib.h>
#include <time.h>
#include <float.h>
#include <math.h>
#include "prng.h"

rand_init_f_t rand_init_f = std_init;
rand_f_t rand_f = std_rand;
rand256i_f_t rand256i_f = std_rand256i;

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

__m256d prng_get_random_double4_in_range(__m256d lower_bound, __m256d upper_bound){
   // TODO: proper SIMD!
   const __m256d rMaxInv = _mm256_set1_pd(1.0/RAND_MAX);
   __m256i ri = rand256i_f();
   __m128i rs = _mm256_castsi256_si128(ri);
   __m256d rr = _mm256_cvtepi32_pd(rs);
   //return rr;
   __m256d r = _mm256_mul_pd(rr,rMaxInv);
   __m256d fac = _mm256_sub_pd(upper_bound, lower_bound);
   __m256d res = _mm256_fmadd_pd(r,fac, lower_bound);
   return res;
   
   // ------------- third try: with bit mangling
   //const __m256i exp = _mm256_set1_epi64x(1023L << 52);
   //const __m256i mask = _mm256_set1_epi64x(0xFFFFFFFF);
   //__m256i r = rand256i_f();
   //r = _mm256_and_si256(r,mask);
   //r = _mm256_slli_epi64(r,21); // 1 lat, 1 tp
   //r = _mm256_or_si256(r,exp); // 1 lat, 2 or 3 throughput
   //__m256d rd = _mm256_castsi256_pd(r);
   //// above: value between 1..2
   //rd = _mm256_sub_pd(rd, _mm256_set1_pd(1));
   //__m256d fac = _mm256_sub_pd(upper_bound, lower_bound);
   //return _mm256_fmadd_pd(rd,fac, lower_bound);

   // -------------- second attempt: still very slow.
   //__m256d rr = _mm256_set_pd(rand_f(), rand_f(), rand_f(), rand_f());
   //const __m256d rMaxInv = _mm256_set1_pd(1.0/RAND_MAX);
   //__m256d r = _mm256_mul_pd(rr,rMaxInv);
   //__m256d fac = _mm256_sub_pd(upper_bound, lower_bound);
   //__m256d res = _mm256_fmadd_pd(r,fac, lower_bound);
   //return res;
   // --------------first attempt: very slow:
   //__m256d r;
   //r[0] = prng_get_random_double_in_range(lower_bound[0],upper_bound[0]);
   //r[1] = prng_get_random_double_in_range(lower_bound[1],upper_bound[1]);
   //r[2] = prng_get_random_double_in_range(lower_bound[2],upper_bound[2]);
   //r[3] = prng_get_random_double_in_range(lower_bound[3],upper_bound[3]);
   //return r;
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
