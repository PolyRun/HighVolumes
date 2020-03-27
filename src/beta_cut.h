#include "poly/volume.h"



/**
 * \brief get initial estimate for the outer (?) ellipsoid around A.
 *
 * \param n dimension of space
 * \param m number of hyperplanes
 * \param A the polytope as specified in volume.h
 * \R2 return value: initial estimate of ellipsoid matrix is R2 * I
 * \Ori return value: initial estimate of central point of ellipsoid
 **/
void initEllipsoid(Polytope *Pol, FT *R2, FT **Ori);
