#ifndef EXAMPLE_POLYTOPES
#define EXAMPLE_POLYTOPES

#include "../../src/volume/volume.h"
#include "../../src/volume/volume_helper.hpp"

struct Solved_Polytope {
    Polytope *polytope;
    FT volume;
};

// Please provide lower_and upper bounds for each dimension of the rectangle
// See function generate_unit_hypercube for an example
struct Solved_Polytope generate_hyperrectangle(int dims, FT *lower_bounds, FT *upper_bounds) {

    int num_constraints = 2*dims;

    Polytope *hyperrectangle = Polytope_new(dims, num_constraints);

    // Lower bounds
    for (int i = 0; i < dims; i++) {
        // Ax >= b means that -Ax <= -b
        Polytope_set_a(hyperrectangle, i, i, -1);
        Polytope_set_b(hyperrectangle, i, -lower_bound[i]);
    }

    // Upper bounds
    for (int i = 0; i < dims; i++) {
        Polytope_set_a(hyperrectangle, i + dims, i, 1);
        Polytope_set_b(hyperrectangle, i, upper_bound[i]);
    }

    FT volume = 1.0;
    for (int i = 0; i < dims; i++) {
        volume *= upper_bound[i] - lower_bound[i];
    }

    struct Solved_Polytope result = { hyperrectangle, volume };
    return result;

}

struct Solved_Polytope generate_unit_hypercube(int dims) {

    FT lower_bounds[dims];
    FT upper_bounds[dims];

    for (int i = 0; i < dims; i++) {
        lower_bounds[i] = -0.5;
        upper_bounds[i] = +0.5;
    }

    return generate_hyperrectangle(dims, lower_bounds, upper_bounds);
    
}

// A cross polytope is the n-dimensional generalisation of a octahedron
// Here it is designed such that its corners are distance 1 away from the origin
// Careful! A cross polytope has 2^n hyperplanes. Don't set n too large!
struct Solved_Polytope generate_cross_polytope(int dims) {

    // A cross-polytope is defined by ± x_1 ± x_2 ... ± x_n <= 1
    int num_constraints = (1 << dims);

    Polytope *cross_polytope = Polytope_new(dims, num_constraints);

    for (int i = 0; i < num_constraints; i++) {
        for (int j = 0; j < dims; j++) {
            // i is a number from 0 to 2^dims - 1
            // We check its j-th bit to decide whether we want x_j or -x_j
            FT plus_minus_1 = 1.0 - 2 * (i & (1 << j));
            Polytope_set_a(cross_polytope, i, j, plus_minus_1);
        }
    }

    for (int i = 0; i < num_constraints; i++) {
        Polytope_set_b(cross_polytope, i, 1);
    }

    // volume of a cross polytope is 2^n / n!
    long n_factorial = 1;
    for (long i = 1; i <= dims; i++) {
        n_factorial *= i;
    }
    FT volume = ((FT) num_constraints) / n_factorial;

    struct Solved_Polytope result = { cross_polytope, volume };
    return result;

}

struct Solved_Polytope generate_solved_polyvest_polytope(int index) {

    Polytope *polytope;

    int error = read_polyvest_p(exp_paths[index], polytope);

    assert(error != 1 && "Cannot generate polyvest polytope. Aborting");

    vol::Polyvest_p reference_polytope(polytope->n, polytope->m);
    polyvest_convert(polytope, &reference_polytope);

    int step_size = 10;

    reference_polytope.Preprocess();
    reference_polytope.EstimateVol(step_size);
    FT volume = (FT) reference_polytope.Volume();

    struct Solved_Polytope result = { polytope, volume };
    return result;

}


#endif