#include <iostream>
#include "../../src/volume/volume_helper.hpp"
#include "test_helpers.hpp"

extern "C" { // must be included C stlye
#include "../../src/volume/preprocess.h"
}


void test_init_against_polyvest(Polytope *P){

    int n = P->n;
    int m = P->m;  
  
    vol::Polyvest_p Q(m, n);
    polyvest_convert(P, &Q);

  
    vec polyvest_ori(n);
    double polyvest_R2;
  
    Q.genInitE(polyvest_R2, polyvest_ori);
  
    FT R2;
    FT *ori;
    init_ellipsoid(P, &R2, &ori);

    // this cast is legit according to
    // https://stackoverflow.com/questions/2923272/how-to-convert-vector-to-array
    FT *polyvest_ori_array = &polyvest_ori[0];
    FT *highvolumes_ori_array = &ori[0];
    FT diff_ori = frobenius(polyvest_ori_array, highvolumes_ori_array, n, 1);
    FT diff_r2 = frobenius(&R2, &polyvest_R2, 1, 1);

    assert(diff_r2 >= 0 && diff_ori >= 0 &&
           "returned ellipsoid does not consist of reals");

    assert(diff_r2 < 0.1 && diff_ori < 0.1 &&
           "difference in frobenius norm is not small enough!");


    std::cout << "PASSED" << std::endl;
    

#ifdef TEST_MSG
    std::cout << "---------------- POLYVEST\n";
    std::cout << "initial R2: " << polyvest_R2 << std::endl;
    std::cout << "initial ori:" << std::endl;
    for (int i = 0; i < n; i++){
        cout << polyvest_ori(i) << " ";
    }
    std::cout << endl;
    
    std::cout << "----------------- Highvolumes\n";
    std::cout << "initial R2: " << R2 << std::endl;
    std::cout << "initial ori:" << std::endl;
    for (int i = 0; i < n; i++){
        cout << ori[i] << " ";
    }
    std::cout << endl;
#endif  
}



int main(){

    
    std::cout << "\n-------------- TEST INIT EXAMPLE POLYTOPES:\n";
    
    Polytope *P;
    for (int i = 0; i < 33; i++){

        std::cout << "TESTING " << exp_paths[i] << std::endl;
        
        int err = read_polyvest_p(exp_paths[i], &P);
        assert(!err);

        test_init_against_polyvest(P);
        
        Polytope_free(P);
    }

    return 0;
    
}


