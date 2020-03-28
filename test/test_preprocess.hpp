
#include "../polyvest/vol.h"
#include <cmath>

extern "C" { // must be included C stlye
#include "../src/poly/volume.h"
}



#ifndef TEST_PREPROCESS_H
#define TEST_PREPROCESS_H


#define EPS 1e-14


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
            assert(std::isfinite(a1) && std::isfinite(a2));
            FT a_diff = a1 - a2;
            a_diff *= a_diff;
            if (a_diff > EPS){
                res.first += a_diff;
            }
        }
        FT b_diff = Polytope_get_b(P, i) - Q->b(i);
        b_diff *= b_diff;
        if (b_diff > EPS){
            res.second += b_diff;
        }
    }

    return res;
    
}


#endif
