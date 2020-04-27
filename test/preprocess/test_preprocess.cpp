#include <iostream>
#include "../../src/volume/volume_helper.hpp"
#include "../test_helpers.hpp"

extern "C" { // must be included C stlye
#include "../../src/volume/volume.h"
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

void test_preprocess_generic_box(const int n, Body_T* type_, void* box, void* box_out) {
   FT det;
   void* body_in[1] = {box};
   void* body_out[1] = {box_out};
   Body_T* type[1] = {type_};

   type_->print(box);
   preprocess_generic(n, 1, (const void**) body_in, (void**) body_out, (const Body_T**) type, &det);
   
   type_->print(box_out);
   
   FT* x = (FT*)aligned_alloc(32, n*sizeof(FT));
   void* cache = aligned_alloc(32, type_->cacheAlloc(box_out));
   for(int i=0;i<n;i++) {x[i]=0;}
   type_->cacheReset(box_out,x,cache);
   // check origin 0 is inside:
   assert(type_->inside(box_out, x));
   
   // Since the box was alligned, we know it is still alligned
   // this is because the cutting planes were all alligned
   // So now we measure from the origin to the sides.
   // then we know how far the sides are out
   // this way, we can find the corner furthest from the origin
   // this cannot be too far out
   FT d2 = 1;// squared distance of furthest point
   for(int i=0;i<n;i++) {
       FT t0,t1;
       type_->intersectCoord(box_out, x, i, &t0, &t1, cache);
       FT tmax = std::max(-t0,t1);
       d2 += tmax*tmax;
       assert(t0 <= -1 && t1 >= 1 && "walls do not cut inner ellipse");
       std::cout << t0 << " " << t1 << "\n";
   }
   std::cout << "d2 " << d2 << " vs " << (4*n*n)<< "\n";
   assert(d2 <= 4.0*n*n && "box not outside ellipse");
   
   free(x);
   free(cache);
   std::cout << "det: " << det << std::endl;
}
   
void test_preprocess_generic() {
    std::cout << "\n ----------- TEST GENERIC PREPROCESSING:\n";
   
    {
        const int n = 10;
	Polytope* box = Polytope_new_box(n,0.5);
	Polytope* box_out = Polytope_new_box(n,1.0);

	test_preprocess_generic_box(n, &Polytope_T, box, box_out);

        Polytope_T.free(box);
        Polytope_T.free(box_out);
    }
    {
        const int n = 10;
        PolytopeT* box = PolytopeT_new_box(n,0.5);
        PolytopeT* box_out = PolytopeT_new_box(n,1.0);

        test_preprocess_generic_box(n, &PolytopeT_T, box, box_out);

        PolytopeT_T.free(box);
        PolytopeT_T.free(box_out);
    }
    
    {
        const int n = 10;
        FT det;
        Ellipsoid* e1 = Ellipsoid_new_with_T(n);
        Ellipsoid* e2 = Ellipsoid_new_with_T(n);
        for(int i=0; i<n; i++) {
            e1->a[i] = (i==1)*10;//prng_get_random_double_in_range(-0.1,0.1);
            e2->a[i] = (i==1)*10;//prng_get_random_double_in_range(-0.1,0.1);
            FT* Ai1 = Ellipsoid_get_Ai(e1,i);
            Ai1[i] = prng_get_random_double_in_range(0.1,0.2);
            FT* Ai2 = Ellipsoid_get_Ai(e2,i);
            Ai2[i] = prng_get_random_double_in_range(0.1,0.2);
        }
        void* body_in[2] = {e1, e2};
        Ellipsoid* e1_out = Ellipsoid_new_with_T(n);
        Ellipsoid* e2_out = Ellipsoid_new_with_T(n);
        void* body_out[2] = {e1_out, e2_out};
        Body_T* type[2] = {&Ellipsoid_T, &Ellipsoid_T};

        preprocess_generic(n, 2, (const void**) body_in, (void**) body_out, (const Body_T**) type, &det);
        
	std::cout << "e1_out:\n";
        Ellipsoid_T.print(e1_out);
	std::cout << "e2_out:\n";
        Ellipsoid_T.print(e2_out);
        
        // check that centers allign:
	for(int i=0;i<n;i++) {
	   assert(e1_out->a[i] == e2_out->a[i]);
	}
        
        // eval test:
        FT* x = (FT*)aligned_alloc(32, n*sizeof(FT));
	for(int i=0;i<n;i++) {
	    {// inner ellipsoid
                for(int j=0;j<n;j++) {x[j]=(i==j);}// init vector
	        FT eval1 = Ellipsoid_eval(e1_out, x);
	        FT eval2 = Ellipsoid_eval(e2_out, x);
		std::cout << "inner: " << eval1 << " " << eval2 << "\n";
		assert(eval1 <= 1.0);
		assert(eval2 <= 1.0);
            }
	    {// outer ellipsoid
                for(int j=0;j<n;j++) {x[j]=(i==j)*2*n;}// init vector
	        FT eval1 = Ellipsoid_eval(e1_out, x);
	        FT eval2 = Ellipsoid_eval(e2_out, x);
		std::cout << "outer: " << eval1 << " " << eval2 << "\n";
		assert(eval1 >= 1.0);
		assert(eval2 >= 1.0);
            }
	}

	std::cout << "det: " << det << std::endl;
        
	free(x);
	Ellipsoid_T.free(e1);
	Ellipsoid_T.free(e2);
	Ellipsoid_T.free(e1_out);
	Ellipsoid_T.free(e2_out);
    }

    {
        const int n = 10;
        FT det;
	Polytope* box = Polytope_new_box(n,1.0);
        Ellipsoid* e = Ellipsoid_new_with_T(n);
        for(int i=0; i<n; i++) {
            e->a[i] = (i==0)*1.7; // set on plane of box
            FT* Ai = Ellipsoid_get_Ai(e,i);
            Ai[i] = 2.0;// make quite small
	    if(i==0) {Ai[i] = 1;}// stretch quite far
        }

        void* body_in[2] = {box, e};
	Polytope* box_out = Polytope_new_box(n,1.0);
        Ellipsoid* e_out = Ellipsoid_new_with_T(n);
        void* body_out[2] = {box_out, e_out};
        Body_T* type[2] = {&Polytope_T, &Ellipsoid_T};

        preprocess_generic(n, 2, (const void**) body_in, (void**) body_out, (const Body_T**) type, &det);
        
	std::cout << "det: " << det << std::endl;

	std::cout << "e_out:\n";
        Ellipsoid_T.print(e_out);
	std::cout << "box_out:\n";
	Polytope_T.print(box_out);

        // inside test:
        FT* x = (FT*)aligned_alloc(32, n*sizeof(FT));
	for(int i=0;i<n;i++) {
	    {// inner ellipsoid
                for(int j=0;j<n;j++) {x[j]=(i==j);}// init vector
                bool i1 = Polytope_T.inside(box_out,x);
                bool i2 = Ellipsoid_T.inside(e_out,x);
		assert(i1 && i2 && "inner ellipsoid points must be in both");
	    }
	    {// outer ellipsoid
                for(int j=0;j<n;j++) {x[j]=(i==j)*2*n;}// init vector
                bool i1 = Polytope_T.inside(box_out,x);
                bool i2 = Ellipsoid_T.inside(e_out,x);
		std::cout << "outer " << i1 << " " << i2 << "\n";
		assert((!i1 && !i2) && "outer ellipsoid points must be outside at least one");
            }
	}

	Polytope_T.free(box);
	Polytope_T.free(box_out);
	Ellipsoid_T.free(e);
	Ellipsoid_T.free(e_out);
    }
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
    std::cout << "\n-------------- TEST PREPROCESS RANDOM POLYTOPES:\n"
              << ntests << " random polytopes of dim " << dim << " with " << nconstraints << " constraints:\n";
    test_preprocess_random_polytopes(ntests, dim, nconstraints);
    
    
    std::cout<< "TESTS COMPLETE.\n";
    
}
