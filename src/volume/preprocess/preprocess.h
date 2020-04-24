#include "../ft.h"
#include "../matrix/matrix.h"
#include "../polytope/polytope.h"
#include "../ellipsoid/ellipsoid.h"
#include "../linalg/linalg.h"



/**
 * \brief get initial estimate for the outer (?) ellipsoid around A.
 *
 * \param n dimension of space
 * \param m number of hyperplanes
 * \param A the polytope as specified in volume.h
 * \R2 return value: initial estimate of ellipsoid matrix is R2 * I
 * \Ori return value: initial estimate of central point of ellipsoid
 **/
void init_ellipsoid(const Polytope *Pol, FT *R2, FT **Ori);





/**
 * \brief transform input polytope s.t. it fits in unit ball and has 1/2n * unit ball inside it
 * \param P: input polytope
 * \param Q: output polytope
 * \param det: output: the determinant of the linear transformation applied to P in order to get Q
 **/
void preprocess(Polytope *P, Polytope **Q, FT *det);


/**
... same as before plus
 * \param c1: counts number of cuts applied
 * \param c2: counts number of iterations in first inner loop (where we check ori)
 * \param c3: counts number of iterations in second inner loop (where we check the inner ellipsoid)
 * \param c4: counts number of times we compute the cut in the first inner loop
 **/
void preprocess_opcount(Polytope *P, FT R2, FT *ori, Polytope **Q, FT *det, int *iterations, int *loopone, int *looptwo, int *breakcond);
