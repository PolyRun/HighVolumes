#ifndef HEADER_DOTPRODUCT_H
#define HEADER_DOTPRODUCT_H

#include "../linalg.h"

// Reference impl:
FT dotProduct_ref(const FT* u, const FT* v, const int n);

// scalar, 2 accumulators
FT dotProduct_2acc(const FT* u, const FT* v, const int n);

// First vectorization:
FT dotProduct_vec1(const FT* u, const FT* v, const int n);
// comment:
//  seems to be slower than ref, takes more than 2x time

FT dotProduct_auto1(const FT* u, const FT* v, const int n);
FT dotProduct_auto2(const FT* u, const FT* v, const int n);

// ------------------------------ squared Norm, i.e. x^T * x
FT squaredNorm_ref(const FT* v, const int n);


#endif // HEADER_DOTPRODUCT_H
