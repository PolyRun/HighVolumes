#include "volume_examples.hpp"

Solved_Body* generate_body(const std::string &generator) {
    if(generator.compare("cube_10") == 0) {
       return generate_unit_hypercube(10);	
    } else if(generator.compare("cube_20") == 0) {
       return generate_unit_hypercube(20);	
    }
}

Solved_Body* generate_hyperrectangle(int dims, FT *lower_bounds, FT *upper_bounds) {
    int num_constraints = 2*dims;

    PolytopeT *hyperrectangle = PolytopeT_new(dims, num_constraints);

    // Lower bounds
    for (int i = 0; i < dims; i++) {
        // Ax >= b means that -Ax <= -b
        PolytopeT_set_a(hyperrectangle, i, i, -1);
        PolytopeT_set_b(hyperrectangle, i, -lower_bounds[i]);
    }

    // Upper bounds
    for (int i = 0; i < dims; i++) {
        PolytopeT_set_a(hyperrectangle, i + dims, i, 1);
        PolytopeT_set_b(hyperrectangle, i, upper_bounds[i]);
    }

    FT volume = 1.0;
    for (int i = 0; i < dims; i++) {
        volume *= upper_bounds[i] - lower_bounds[i];
    }

    Solved_Body* result = new Solved_Body(1);
    result->body[0] = hyperrectangle;
    result->type[0] = &PolytopeT_T;
    result->volume = volume;
    return result;
}

Solved_Body* generate_unit_hypercube(int dims) {
    FT lower_bounds[dims];
    FT upper_bounds[dims];

    for (int i = 0; i < dims; i++) {
        lower_bounds[i] = -0.5;
        upper_bounds[i] = +0.5;
    }

    return generate_hyperrectangle(dims, lower_bounds, upper_bounds);
}

Solved_Body* generate_cross_polytope(int dims) {

    // This cross-polytope is defined by ± x_1 ± x_2 ... ± x_n <= 1
    int num_constraints = (1 << dims);

    PolytopeT *cross_polytope = PolytopeT_new(dims, num_constraints);

    for (int i = 0; i < num_constraints; i++) {
        for (int j = 0; j < dims; j++) {
            // i is a number from 0 to 2^dims - 1
            // We check its j-th bit to decide whether we want x_j or -x_j
            FT plus_minus_1 = 1.0 - 2 * ((i & (1 << j)) >> j);
            PolytopeT_set_a(cross_polytope, i, j, plus_minus_1);
        }
    }

    for (int i = 0; i < num_constraints; i++) {
        PolytopeT_set_b(cross_polytope, i, 1);
    }

    // Volume of a cross polytope is 2^n / n!
    long n_factorial = 1;
    for (long i = 1; i <= dims; i++) {
        n_factorial *= i;
    }
    FT volume = ((FT) num_constraints) / n_factorial;

    int bcount = 1;
    Solved_Body *result = new Solved_Body(bcount);
    result->body[0] = cross_polytope;
    result->type[0] = &PolytopeT_T;
    result->volume = volume;
    return result;

}

Solved_Body* generate_simplex(int dims) {

    // We define the simplex as x_i >= 0 for all i
    // and one more constraint x_1 + ... + x_n <= 2
    // If the constraint was "<= 1", then its volume would be 1 / n!
    // But by scaling every side by two, its volume becomes 2^n / n!
    // This prevents the volume from going to zero too fast
    int num_constraints = dims + 1;

    PolytopeT *simplex = PolytopeT_new(dims, num_constraints);

    // x_i >= 0
    for (int i = 0; i < dims; i++) {
        PolytopeT_set_a(simplex, i, i, -1);
        PolytopeT_set_b(simplex, i, 0);
    }

    // x_1 + ... x_n <= 2
    for (int i = 0; i < dims; i++) {
        PolytopeT_set_a(simplex, dims, i, 1);
    }
    PolytopeT_set_b(simplex, dims, 2);

    // Again, volume = 2^n / n!, just like the cross polytope!
    long n_factorial = 1;
    for (long i = 1; i <= dims; i++) {
        n_factorial *= i;
    }
    FT volume = ((FT) num_constraints) / n_factorial;

    int bcount = 1;
    Solved_Body *result = new Solved_Body(bcount);
    result->body[0] = simplex;
    result->type[0] = &PolytopeT_T;
    result->volume = volume;
    return result;

}

Solved_Body* generate_ellipsoid(int dims, FT *lower_bounds, FT *upper_bounds) {

    Ellipsoid *ellipsoid = Ellipsoid_new(dims);

    FT volume = 1.0;

    for (int i = 0; i < dims; i++) {
        FT radius_i = (upper_bounds[i] - lower_bounds[i]) / 2;
        FT midpoint = lower_bounds[i] + radius_i;

        ellipsoid->a[i] = midpoint;
        
        // Instead of x_1^2 + ... + x_n^2 <= 1
        // we want that (x_1/r_1)^2 + ... + (x_n/r_n)^2 <= 1
        FT lambda_i = 1 / (radius_i * radius_i);
        Ellipsoid_set_a(ellipsoid, i, i, lambda_i);

        // Over all i's we calculate det(T) here, where T = A^{-1}
        volume *= radius_i * radius_i;
    }

    int bcount = 1;
    Solved_Body *result = new Solved_Body(bcount);
    result->body[0] = ellipsoid;
    result->type[0] = &Ellipsoid_T;
    result->volume = volume;
    return result;

}

Solved_Body* generate_unit_ball(int dims) {

    FT lower_bounds[dims];
    FT upper_bounds[dims];

    for (int i = 0; i < dims; i++) {
        lower_bounds[i] = -0.5;
        upper_bounds[i] = +0.5;
    }

    return generate_ellipsoid(dims, lower_bounds, upper_bounds);

}
//
//// This function is out of order, because read_polyvest_p() doesn't exist anymore
//// We need to reimplement it anyway to return PolytopeT instead of Polytope
///*
//struct Solved_Body generate_solved_polyvest_polytope(int index) {
//
//    PolytopeT *polytope;
//
//    int error = read_polyvest_p(exp_paths[index], polytope);
//
//    assert(error != 1 && "Cannot generate polyvest polytope. Aborting");
//
//    vol::Polyvest_p reference_polytope(polytope->n, polytope->m);
//    polyvest_convert(polytope, &reference_polytope);
//
//    int step_size = 10;
//
//    reference_polytope.Preprocess();
//    reference_polytope.EstimateVol(step_size);
//    FT volume = (FT) reference_polytope.Volume();
//
//    int bcount = 1;
//    Body_T *body[bcount] = { polytope };
//    struct Solved_Body result = { body, bcount, volume };
//    return result;
//
//}
//*/ 

