#include <iostream>

#include "../test_helpers.hpp"
#include "../../src/volume/volume_helper.hpp"
#include "../../src/random/prng.h"


void test_matrix_invert_pdsym(){

    Matrix *M = Matrix_new(3, 3);
    Matrix *Minv = Matrix_new(3, 3);

    // make M a symmetric, positive-definite matrix
    FT *M1 = Matrix_get_row(M, 0);
    M1[0] = 2; M1[1] = -1; M1[2] = 0;
    FT *M2 = Matrix_get_row(M, 1);
    M2[0] = -1; M2[1] = 2; M2[2] = -1;
    FT *M3 = Matrix_get_row(M, 2);
    M3[0] = 0; M3[1] = -1; M3[2] = 2;

    Matrix_invert_pdsym(M, Minv);

    
    for (int i = 0; i < 3; i++){
        FT *Minvi = Matrix_get_row(Minv, i);
        for (int j = 0; j < 3; j++){
            FT *Mj = Matrix_get_row(M, j);
            FT res = dotProduct(Minvi, Mj, 3);
            //cout << res << "\t";
            if (i == j){
                assert(abs(res - 1) < EPS
                       && "M*Minv[i][i] = 1");
            }
            else {
                assert(abs(res) < EPS
                       && "M*Minv[i][j] = 0");
            
            }
        }
        //cout << "\n";
    }

    
}



void test_matrix_invert_pdsym_random(int n){

    // generate random L matrix
    Matrix *L = Matrix_new(n, n);
    
    for (int i = 0; i < n; i++){
        for (int j = 0; j <= i; j++){
            Matrix_set(L, i, j, prng_get_random_double_in_range(2, 20));
        }
        for (int j = i+1; j < n; j++){
            Matrix_set(L, i, j, 0);
        }
    }

    // compute Ellipsoid's T matrix as L*LT
    Matrix *M = Matrix_new(n,n);

    for (int i = 0; i < n; i++){
        FT *Mi = Matrix_get_row(M, i);
        FT *Li = Matrix_get_row(L, i);
        for (int j = 0; j < n; j++){
            FT *Lj = Matrix_get_row(L, j);
            Mi[j] = dotProduct(Li, Lj, n);
        }
    }

    // compute the inverse
    Matrix *Minv = Matrix_new(n,n);
    Matrix_invert_pdsym(M, Minv);


    FT maxdif = 0;
    for (int i = 0; i < n; i++){
        FT *Minvi = Matrix_get_row(Minv, i);
        for (int j = 0; j < n; j++){
            FT *Mj = Matrix_get_row(M, j);
            FT res = dotProduct(Minvi, Mj, n);

            assert(std::isfinite(res) && "Matrix_invert_pdsym returned nans\n");

            /*
            if (!std::isfinite(res)){
                Matrix_print((void *) Minv);
                cout << "\n";
                Matrix_print((void *) M);
                cout << "\n";
                for (int k = 0; k < n; k++){
                    cout << Matrix_get(L, k, k) << "\n";
                }
                return false;
            }
            */            
            maxdif = (i == j) ?
                maxdif = max(maxdif, abs(res - 1)) :
                abs(res);
        }
    }

    cout << "maxdif: " << maxdif << "\n";

    assert(maxdif < 1e-5 && "M * Minv = I\n");
    
}






int main(){
    
    std::cout << "\n-------------- TEST MATRIX_INVERT_PDSYM \n";

    test_matrix_invert_pdsym();
   
    prng_init();
    
    for (int i = 20; i < 50; i++) {
        test_matrix_invert_pdsym_random(i);
    }


    std::cout<< "TESTS COMPLETE.\n";
    
}
