#ifndef EXAMPLE_POLYTOPES
#define EXAMPLE_POLYTOPES

#include "volume_helper.hpp"

class Solved_Body {
public:
    Solved_Body(const int bcount) : bcount(bcount) {
       body = (void**)malloc(bcount*sizeof(void*));
       type = (Body_T**)malloc(bcount*sizeof(Body_T*));
    }
    ~Solved_Body() {
       delete body;
       delete type;
    }
    void **body;
    Body_T **type;
    int bcount;
    FT volume;
};

// Please provide lower_and upper bounds for each dimension of the rectangle
// See function generate_unit_hypercube for an example
Solved_Body* generate_hyperrectangle(int dims, FT *lower_bounds, FT *upper_bounds);

// A convenience function for generate_hyperrectangle()
// where lower and upper bounds are -0.5 and 0.5 resp.
Solved_Body* generate_unit_hypercube(int dims);

// A cross polytope is the n-dimensional generalisation of a octahedron
// Here it is designed such that its corners are distance 1 away from the origin
// Careful! A cross polytope has 2^n hyperplanes. Don't set n too large!
Solved_Body* generate_cross_polytope(int dims);

// Here we have the n-dimensional generalisation of a right triangle,
// meaning that all but one side are aligned with the axes
// We therefore have n + 1 constraints
Solved_Body* generate_simplex(int dims);

// Generates an ellipse with predefined ranges for each axis
Solved_Body* generate_ellipsoid(int dims, FT *lower_bounds, FT *upper_bounds);

// A convience function for generate_ellipsoid()
// where lower and upper bounds are -0.5 and 0.5 resp.
Solved_Body* generate_unit_ball(int dims);
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


#endif
