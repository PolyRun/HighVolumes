#ifndef EXAMPLE_POLYTOPES
#define EXAMPLE_POLYTOPES

#include "volume_helper.hpp"

class Solved_Body {
public:
    Solved_Body(const int bcount, const int n) : bcount(bcount), n(n) {
       body = (void**)malloc(bcount*sizeof(void*));
       type = (Body_T**)malloc(bcount*sizeof(Body_T*));
       for(int b=0; b<bcount; b++) {
           body[b] = NULL;
           type[b] = NULL;
       }
    }
    ~Solved_Body() {
       for(int b=0; b<bcount; b++) {
          if(body[b]) {
	     type[b]->free(body[b]);
	  }
       }
       delete body;
       delete type;
    }
    void print() {
       std::cout << "Solved_Body: n="<<n<<", bcount="<<bcount<<", volume="<<volume<<"\n";
       for(int b=0; b<bcount; b++) {
          type[b]->print(body[b]);
       }
    }
    
    Solved_Body* clone();
    
    // x = (L * y + a)*beta
    // smaller beta, body grows
    // det = det(L)
    Solved_Body* transform(const Matrix* L, const FT det, const FT* a, const FT beta);
    Solved_Body* scale(const FT beta);
    Solved_Body* rotate();
    
    void **body;
    Body_T **type;
    int bcount;
    FT volume;
    int n;
    bool is_normalized = false;
};

// ----------------------------------------------------------- Solved_Body Generator
class Solved_Body_Generator {
public:
    Solved_Body_Generator(); // here all the generators are registered!
    Solved_Body* get(std::string name) {
        auto it = generators.find(name);
	if(it!=generators.end()) {
	    return it->second();
	} else {
	    return NULL;
	}
    }
    void add(const std::string &name, const std::string &desc, std::function<Solved_Body*()> gen) {
        generators[name] = gen;
	gen_map_[name] = {name,desc};
    }
    const std::map<std::string, std::pair<std::string, std::string>>& gen_map() {return gen_map_;}
private:
    std::map<std::string, std::function<Solved_Body*()>> generators;
    std::map<std::string, std::pair<std::string,std::string>> gen_map_;
    // map name to {name, desc}
    // used for CLI
};

Solved_Body_Generator* solved_body_generator();


// ----------------------------------------------------------- Generator Functions

// Please provide lower_and upper bounds for each dimension of the rectangle
// See function generate_centered_hypercube for an example
Solved_Body* generate_hyperrectangle(int dims, FT *lower_bounds, FT *upper_bounds);

// A convenience function for generate_hyperrectangle()
// where lower and upper bounds are -r and r resp.
Solved_Body* generate_centered_hypercube(int dims,FT r);

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
// where lower and upper bounds are -r and r resp.
Solved_Body* generate_centered_ball(int dims, FT r);
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
