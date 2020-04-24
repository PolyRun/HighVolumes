#include <iostream>
#include <cassert>

#include "../../src/util/cli.hpp"
#include "../../src/volume/preprocess.h"
#include "../../src/volume/volume_helper.hpp"
#include "../preprocess/test_helpers.hpp" // stealing manuels polytopes


// This test compares the volumes of polyvest with the volumes we calculate

void test_polyvest_example(int index) {

    assert(("Polyvest example index not in range",
            0 <= index && index < NEXAMPLE_POLYTOPES));
    
    Polytope *P;

    int error = read_polyvest_p(exp_paths[index], P);
    if (error == 1) {
        std::cout << "Failed to read polytope " << exp_paths[index] << std::endl;
        return;
    }

    int dims = P->n;
    FT *scaling_factor;
    Polytope *P_normalized;
    
    int nums = 1;
    Polytope *bodies[nums]            = {P};
    Polytope *normalized_bodies[nums] = {P_normalized};
    Body_T   *body_types[nums]        = {&Polytope_T};

    preprocess_ref(dims, nums, (const void **) bodies,
                               (const void **) normalized_bodies,
                               body_types,
                               scaling_factor);
    
    FT inner_radius = 0.00000000001;
    FT outer_radius = 1000000000000;

    FT volume;
    volume = volume_ref(dims, inner_radius, outer_radius, normalized_bodies, body_types);

    std::cout << "Volume of " << exp_paths[index] << ": " << volume << std::endl;

    Polytope_free(P);

}

void test_all_polyvest_examples() {

    for (int i = 0; i < NEXAMPLE_POLYTOPES; i++) {
        test_polyvest_example(i);
    }

}

int main(int argc, char **argv) {
    CLI cli(argc, argv, "End to end test");

    
}