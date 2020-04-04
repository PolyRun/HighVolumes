#include "volume.h"



/**
 * \brief get initial estimate for the outer (?) ellipsoid around A.
 *
 * \param n dimension of space
 * \param m number of hyperplanes
 * \param A the polytope as specified in volume.h
 * \R2 return value: initial estimate of ellipsoid matrix is R2 * I
 * \Ori return value: initial estimate of central point of ellipsoid
 **/
void init_ellipsoid(Polytope *Pol, FT *R2, FT **Ori);




/**
 * \brief transform input polytope s.t. it fits in unit ball and has 1/2n * unit ball inside it
 * \param P: input polytope
 * \param Q: output polytope
 * \param det: output: the determinant of the linear transformation applied to P in order to get Q
 **/
void preprocess(Polytope *P, Polytope **Q, double *det);
