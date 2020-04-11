#include <iostream>
#include "../../src/volume/volume_helper.hpp"
#include "test_helpers.hpp"

extern "C" { // must be included C stlye
#include "../../src/volume/preprocess.h"
}



void test_preprocess_against_polyvest(Polytope *P){

    int n = P->n;
    int m = P->m;  
  
    vol::Polyvest_p Q(m, n);
    polyvest_convert(P, &Q);

    //Q.A.print();
    //Q.b.print();


    Q.Preprocess();

    Polytope *R;
    FT det;
    preprocess(P, &R, &det);


    std::pair<FT, FT> diff = matrix_diff(R, &Q);

    assert(diff.second >= 0 && diff.first >= 0 &&
           "preprocess returned non-real polytope");
    
    assert(diff.first < 0.1 && diff.second < 0.1 &&
           "frobenius norm of preprocessed polytopes is too different");

    std::cout << "PASSED" << std::endl;
        

#ifdef TEST_MSG    
    std::cout << "2-Frobenius of A_P - A_Q: "
              << diff.first << std::endl
              << "2-norm of b_P - b_Q:      "
              << diff.second << std::endl;    
#endif
}

/**
 *\brief test if the scaled polytope contains scaled unit ball B(0, 1/(2n))
 **/
void test_preprocess_circumscription(Polytope *P){
    
    Polytope *R;
    FT det;
    preprocess(P, &R, &det);

    bool is_correct = polytope_contains_scaled_ball(R);

    assert(is_correct &&
           "returned polytope does not include scaled 1-ball!");

    std::cout << "PASSED" << std::endl;
        

    
#ifdef TEST_MSG
    if (is_correct){
        std::cout << "Polytope contains B(0, 1/(2n))" << std::endl;
    }
    else {
        std::cout << "ERROR: Polytope doesn't contain B(0, 1/(2n))" << std::endl;
    }
#endif
    
}


void test_preprocess_example_polytopes(){
    Polytope *P;
    for (int i = 0; i < NEXAMPLE_POLYTOPES; i++){

        std::cout << "TESTING " << exp_paths[i] << std::endl;
        
        int err = read_polyvest_p(exp_paths[i], &P);
        assert(!err);

        test_preprocess_against_polyvest(P);
        test_preprocess_circumscription(P);

    }
    
    Polytope_free(P);
}



// TODO: choosing ntests > 1 for the moment doesn't make sense as my rng seed changes too slow...
void test_preprocess_random_polytopes(int ntests, int dim, int nconstraints){

    // unit ball
    std::vector<double> ell(dim, 1.0);
    Polytope *P;
    
    for (int i = 0; i < ntests; i++){
        
        std::cout << "TESTING RANDOM POLYTOPE " << i << std::endl;

        make_random_poly(ell, nconstraints, &P);

        test_preprocess_against_polyvest(P);
        test_preprocess_circumscription(P);

    }

    Polytope_free(P);
    
}


int main(){    
    std::cout << "\n-------------- TEST PREPROCESS EXAMPLE POLYTOPES:\n";
    test_preprocess_example_polytopes();

    int ntests = 1;
    int dim = 20;
    int nconstraints = 100;
    std::cout << "\n-------------- TEST PREPROCESS EXAMPLE POLYTOPES\n"
              << ntests << " random polytopes of dim " << dim << " with " << nconstraints << " constraints:\n";
    test_preprocess_random_polytopes(ntests, dim, nconstraints);
    
    
}
