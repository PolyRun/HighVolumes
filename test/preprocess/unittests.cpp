#include "test_helpers.hpp"
#include <iostream>
#include "../../src/volume/volume_helper.hpp"



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

    cholesky_ellipsoid(M, L);

    
    Ellipsoid_print(M);
    Matrix_print(L);
    
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

    std::cout<< "TESTS COMPLETE.\n";
    
}
