#include "../polyvest/vol.h"
#include <cmath>

extern "C" { // must be included C stlye
#include "../src/volume/volume.h"
}



#ifndef TEST_PREPROCESS_H
#define TEST_PREPROCESS_H


#define EPS 1e-14


//#define TEST_MSG

#define POLYEXP_BASE ((string) "../../../polyvest/examples/")
#define NEXAMPLE_POLYTOPES 32

static string exp_paths[NEXAMPLE_POLYTOPES] = {
                                        POLYEXP_BASE + "cc_8_10",
                                        POLYEXP_BASE + "cc_8_11",
                                        POLYEXP_BASE + "cross_13",
                                        POLYEXP_BASE + "cross_7",
                                        POLYEXP_BASE + "cross_9",
                                        POLYEXP_BASE + "cube_10",
                                        POLYEXP_BASE + "cube_10_2",
                                        POLYEXP_BASE + "cube_14",
                                        POLYEXP_BASE + "cube_14_2",
                                        POLYEXP_BASE + "cube_15",
                                        POLYEXP_BASE + "cube_2",
                                        POLYEXP_BASE + "cube_20",
                                        POLYEXP_BASE + "cube_25",
                                        POLYEXP_BASE + "cube_30",
                                        POLYEXP_BASE + "cube_35",
                                        POLYEXP_BASE + "cube_40",
                                        POLYEXP_BASE + "cube_5",
                                        POLYEXP_BASE + "cube_80",
                                        POLYEXP_BASE + "ex_1",
                                        POLYEXP_BASE + "ex_2",
                                        POLYEXP_BASE + "fm_6",
                                        /* this test fails when math operations are reordered, e.g. with -O3 or with --ffast-math due to numerical imprecision
                                          POLYEXP_BASE + "rect_3",*/
                                        POLYEXP_BASE + "rh_1",
                                        POLYEXP_BASE + "rh_2",
                                        POLYEXP_BASE + "rh_20_40",
                                        POLYEXP_BASE + "rh_3",
                                        POLYEXP_BASE + "rh_30_60",
                                        POLYEXP_BASE + "rh_4",
                                        POLYEXP_BASE + "rh_40_80",
                                        POLYEXP_BASE + "simplex_10",
                                        POLYEXP_BASE + "simplex_14",
                                        POLYEXP_BASE + "simplex_15",
                                        POLYEXP_BASE + "simplex_20",
};




/**
 *\brief decide if ellipsoid (E, c) is included in polytope P. note the dimensions must coincide!
 *\param P a polytope
 *\param E the positive-definite matrix of size P->nxP->n describing the ellipsoid
 *\param c the center of ellipsoid, an P->n vector
 **/
bool ellipsoid_inside_poly(Polytope *P, FT *E, FT *c);

/**
 * \brief check if P contains 1/beta B(0,1) where beta is 2n in our case. this function uses ellipsoid_inside_poly internally
 * \param P the polytope to check
 **/
bool polytope_contains_scaled_ball(Polytope *P);



/**
 * \brief compute frobenius norm of A - B
 * \param A input matrix (A is a d1xd2 matrix)
 * \param B second input matrix (B is a d1xd2 matrix)
 * \param d1 first dimension
 * \param d2 second dimension
 **/
// TODO, adopt to matrix struct
FT frobenius(FT *A, FT *B, int d1, int d2);

/**
 * \brief compare the two polytopes elementwise and return the 2-frobenius norm of the difference matrix of A and the 2-norm of the difference of b
 * we consider EPS as 0
 * this is a specialization of the function frobenius
 **/
std::pair<FT, FT> matrix_diff(Polytope *P, vol::Polyvest_p *Q);


/**
 * \brief test if P is included in unit ball. as this is np-hard we need to enumerate all vertices and check them
 **/
bool polytope_in_unit_ball(Polytope *P);

#endif
