#include "../test_helpers.hpp"
#include <iostream>
#include <iomanip>
#include "../../src/volume/volume_helper.hpp"
#include "../../src/random/prng.h"


void test_cholesky(int n){

    Ellipsoid *M = Ellipsoid_new_with_T(3);
    Matrix *L = Matrix_new(3, 3);
    memset(L->data, 0, 3*3*sizeof(FT));

    // make M a symmetric, positive-definite matrix
    FT *M1 = Ellipsoid_get_Ti(M, 0);
    M1[0] = 2; M1[1] = -1; M1[2] = 0;
    FT *M2 = Ellipsoid_get_Ti(M, 1);
    M2[0] = -1; M2[1] = 2; M2[2] = -1;
    FT *M3 = Ellipsoid_get_Ti(M, 2);
    M3[0] = 0; M3[1] = -1; M3[2] = 2;

    int err = cholesky_ellipsoid(M, L);
    assert(!err && "no cholesky err");

    
    for (int i = 0; i < 3; i++){
        FT *Lrowi = Matrix_get_row(L, i);
        for (int j = 0; j < 3; j++){
            // row of L is col of LT
            FT *LTcolj = Matrix_get_row(L, j);
            assert(abs(
                       dotProduct(Lrowi, LTcolj, 3) -
                       Ellipsoid_get_Ti(M, i)[j]) < EPS
                   && "L * LT != M\n");
        }
    }

}




void test_cholesky_random(int n){

    // generate random L matrix
    Matrix *L = Matrix_new(n, n);

    for (int i = 0; i < n; i++){
        // a larger range incurs larger relative errors...
        // also cholesky can fail if range gets too large
        for (int j = 0; j <= i; j++){
            Matrix_set(L, i, j, prng_get_random_double_in_range(2, 20));
        }
        for (int j = i+1; j < n; j++){
            Matrix_set(L, i, j, 0);
        }
    }

    // compute Ellipsoid's T matrix as L*LT
    Ellipsoid *E = Ellipsoid_new_with_T(n);

    for (int i = 0; i < n; i++){
        FT *Ti = Ellipsoid_get_Ti(E, i);
        FT *Li = Matrix_get_row(L, i);
        for (int j = 0; j < n; j++){
            FT *Lj = Matrix_get_row(L, j);
            Ti[j] = dotProduct(Li, Lj, n);
        }
    }

    // compute the cholesky decomposition
    Matrix *R = Matrix_new(n,n);
    int err = cholesky_ellipsoid(E, R);

    assert(err == 0 && "cholesky failed due to numerical problems");
    /*
    if (err){
        cout << "cholesky failed due to numerical problems!\n"
             << "n: " << n << "\n";
        Ellipsoid_print(E);
        Matrix_print(L);
    }
    */
    FT maxdif = 0;
    for (int i = 0; i < n; i++){
        FT *Rrowi = Matrix_get_row(R, i);
        for (int j = 0; j < n; j++){
            // row of L is col of LT
            FT *RTcolj = Matrix_get_row(R, j);
            FT val = dotProduct(Rrowi, RTcolj, n);
            if(!std::isfinite(val) && "nannannan"){
                Matrix_print(R);
                cout << "error is: " << err << "\n";
                assert(false);
            }
            maxdif = max(maxdif,
                         abs(
                             val -
                             Ellipsoid_get_Ti(E, i)[j]));
            
        }
    }

    /*
    cout << "maxdif " << maxdif << "\n"
         << "n " << n << "\n";
    */
    assert(maxdif < 1e-10 && "R * RT != E.T\n");
                        

}





int main(){
    
    std::cout << "\n-------------- TEST IN BALL\n";

    /*
    for (int dim = 3; dim < 12; dim++){
        double rad = 1.0/sqrt((double) dim) - EPS;
        double rad2 = rad + 0.01;
        std::cout << rad << std::endl << rad2 << std::endl;
    
        Polytope *P = Polytope_new_box(dim, rad);
        //std::cout << P << std::endl;
        Polytope *Q = Polytope_new_box(dim, rad2);
        //std::cout << Q << std::endl;
    
        assert(polytope_in_unit_ball(P));
        assert(!polytope_in_unit_ball(Q));
    }
    */

    test_cholesky(3);

    prng_init();
    
    for (int i = 20; i < 50; i++){
        test_cholesky_random(i);
    }
    

    std::cout<< "TESTS COMPLETE.\n";
    
}
