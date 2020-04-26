#include <iostream>
#include "../../src/volume/volume_helper.hpp"
#include "../test_helpers.hpp"

extern "C" { // must be included C stlye
#include "../../src/volume/volume.h"
}


std::string path_from_exec = "";

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

    FT *polyvest_ori_array = polyvest_ori.memptr();
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

    
    std::cout << "\n-------------- TEST INIT EXAMPLE POLYTOPES:\n";
    
    Polytope *P;
    for (int i = 0; i < NEXAMPLE_POLYTOPES; i++){

        std::cout << "TESTING " << path_from_exec + exp_paths[i] << std::endl;
        
        int err = read_polyvest_p(path_from_exec + exp_paths[i], &P);
        assert(!err);

        test_init_against_polyvest(P);
        
        Polytope_free(P);
    }

    
    std::cout<< "TESTS COMPLETE.\n";

    
}


