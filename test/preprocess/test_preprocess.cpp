#include <iostream>
#include "../../src/volume/volume_helper.hpp"
#include "test_helpers.hpp"

extern "C" { // must be included C stlye
#include "../../src/volume/preprocess.h"
}

std::string path_from_exec = "";

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

        std::cout << "TESTING " << path_from_exec + exp_paths[i] << std::endl;
        
        int err = read_polyvest_p(path_from_exec + exp_paths[i], &P);
        assert(!err &&
               "couldn't read example polytope");

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

void test_preprocess_generic() {
    std::cout << "\n ----------- TEST GENERIC PREPROCESSING:\n";
   
    {
        const int n = 20;
        FT det;
	Polytope* box = Polytope_new_box(n,1.0);
        void* body_in[1] = {box};
	Polytope* box_out = Polytope_new_box(n,1.0);
        void* body_out[1] = {box_out};
        Body_T* type[1] = {&Polytope_T};

        preprocess_ref(n, 1, (const void**) body_in, (void**) body_out, (const Body_T**) type, &det);

	std::cout << "det: " << det << std::endl;
    }
    
    {
        const int n = 20;
        FT det;
	Polytope* box = Polytope_new_box(n,1.0);
        Ellipsoid* e = Ellipsoid_new(n);
        for(int i=0; i<n; i++) {
            e->a[i] = prng_get_random_double_in_range(-0.1,0.1);
            FT* Ai = Ellipsoid_get_Ai(e,i);
            Ai[i] = prng_get_random_double_in_range(1.0,2.0);
        }
        void* body_in[2] = {box, e};
	Polytope* box_out = Polytope_new_box(n,1.0);
        Ellipsoid* e_out = Ellipsoid_new(n);
        void* body_out[2] = {box_out, e_out};
        Body_T* type[2] = {&Polytope_T, &Ellipsoid_T};

        preprocess_ref(n, 2, (const void**) body_in, (void**) body_out, (const Body_T**) type, &det);

	std::cout << "det: " << det << std::endl;
    }
    assert(false && "done.");
}

int main(int argc, char **argv){
    CLI cli(argc, argv, "test preprocess");
    CLIFunctionsVolume cliFun(cli);
    cliFun.preParse();
    if (!cli.parse()){
        return -1;
    }
    cliFun.postParse();

    std::string path = cli.getPath();
    reverse(path.begin(), path.end());
    size_t pos = path.find('/');
    // the executable is not in the current directory
    if(pos != std::string::npos){
        reverse(path.begin(), path.end());
        path_from_exec = path.substr(0, path.length() - pos);
    }
            
    test_preprocess_generic();

    std::cout << "\n-------------- TEST PREPROCESS EXAMPLE POLYTOPES:\n";
    test_preprocess_example_polytopes();

    int ntests = 1;
    int dim = 20;
    int nconstraints = 100;
    std::cout << "\n-------------- TEST PREPROCESS EXAMPLE POLYTOPES\n"
              << ntests << " random polytopes of dim " << dim << " with " << nconstraints << " constraints:\n";
    test_preprocess_random_polytopes(ntests, dim, nconstraints);
    
    
    std::cout<< "TESTS COMPLETE.\n";
    
}
