#include "test_helpers.hpp"
#include <iostream>
#include "../../src/volume/volume_helper.hpp"

int main(){

    std::cout << "\n-------------- TEST IN BALL\n";

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

    std::cout<< "TESTS COMPLETE.\n";
    
}
