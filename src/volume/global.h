#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <immintrin.h>
#include <assert.h>

#include "../random/prng.h"
#include "cholesky.h"



typedef double FT;
#define FT_EPS DBL_EPSILON
#define FT_MAX DBL_MAX
#define FT_MIN DBL_MIN



// dotProduct:
//   assume no memory allignment!
typedef FT (*dotProduct_f_t)(const FT* u, const FT* v, const int n);
extern dotProduct_f_t dotProduct;
#include "dotProduct/dotProduct.h"

// intersect line x + t*d
//           with ball(0,r)
// assume x is in ball
// return t0,t1 for intersections
void Ball_intersect(const int n, const FT r, const FT* x, const FT* d, FT* t0, FT* t1);

void Ball_intersectCoord(const int n, const FT r, const FT* x, const int d, FT* t0, FT* t1);

    
// calculate volume exactly for n-dim ball with radius r
FT Ball_volume(const int n, const FT r);

// given n elements of b bytes, want to get smallest multiple of 32 bytes that fits this.
int ceil_cache(const int n, const int b);

