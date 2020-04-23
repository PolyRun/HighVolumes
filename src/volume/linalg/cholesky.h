
#ifndef CHOLESKY_H
#define CHOLESKY_H

#include "../volume.h"

/**
 * \brief cholesky factorization
 * \param A a symmetric, positive definite matrix
 * \param Trans will be lower diagonal matrix holding L
 * \param n the number of rows and cols of A
 **/
int cholesky(FT *A, FT *Trans, int n);
int cholesky_matrix(const Matrix *M, Matrix *L);
int cholesky_ellipsoid(const Ellipsoid *T, Matrix *L);


#endif
