
#include "../polyvest/vol.h"
#include <cmath>

extern "C" { // must be included C stlye
#include "../src/poly/volume.h"
}



#ifndef TEST_PREPROCESS_H
#define TEST_PREPROCESS_H


#define EPS 1e-14



/**
 *\brief decide if ellipsoid (E, c) is included in polytope P. note the dimensions must coincide!
 *\param P a polytope
 *\param E the positive-definite matrix of size P->nxP->n describing the ellipsoid
 *\param c the center of ellipsoid, an P->n vector
 **/
bool ellipsoid_inside_poly(Polytope *P, FT *E, FT *c);


/**
 * \brief compare the two polytopes elementwise and return the 2-frobenius norm of the difference matrix of A and the 2-norm of the difference of b
 * we consider EPS as 0
 **/
std::pair<FT, FT> matrix_diff(Polytope *P, vol::Polyvest_p *Q){

    int n = P->n;
    int m = P->m;

    std::pair<FT, FT> res = make_pair(0.0, 0.0);
    
    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            FT a1 = Polytope_get_a(P, i, j);
            FT a2 = Q->A(i, j);
            FT a_diff = a1 - a2;
            a_diff *= a_diff;
            if (a_diff > EPS){
                res.first += a_diff;
            }
            if (!std::isfinite(a1) || !std::isfinite(a2)){
                res = std::make_pair(-1, -1);
                break;
            }
        }
        FT b1 = Polytope_get_b(P, i);
        FT b2 = Q->b(i);
        FT b_diff = b1 - b2; 
        b_diff *= b_diff;
        if (b_diff > EPS){
            res.second += b_diff;
        }        
        if (!std::isfinite(b1) || !std::isfinite(b2)){
            res = std::make_pair(-1, -1);
            break;
        }
    }

    return res;
    
}


#endif
