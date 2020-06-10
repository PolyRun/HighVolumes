#include "volume_examples.hpp"
#include <cstdlib>


Solved_Body_Generator::Solved_Body_Generator() {
    // cube
    std::vector<int> cube_n = {3,10,20,40,60,80,100};
    for(int n : cube_n) {
       std::string nstr = std::to_string(n);
       add("cube_r1.0_"+nstr, "basic "+nstr+"-dim cube, centered, side=2 [normalized]", [n]() {
           Solved_Body* sb = generate_centered_hypercube(n,1.0);
           sb->is_normalized = true;
           return sb;
       });
       add("cube_rot_r1.0_"+nstr, "basic "+nstr+"-dim cube, centered, randomly rotated, side=2 [normalized]", [n]() {
           Solved_Body* s = generate_centered_hypercube(n,1.0);
           Solved_Body* sb = s->rotate();
	   delete s;
	   sb->is_normalized = true;
           return sb;
       });
       add("cube_r1.2_"+nstr, "basic "+nstr+"-dim cube, centered, side=2.4 [normalized]", [n]() {
           Solved_Body* sb = generate_centered_hypercube(n,1.2);
           sb->is_normalized = true;
           return sb;
       });
    }
    
    // cross_polytope
    std::vector<int> cross_n = {3,4,5,6,7,8,9,10,11,12,13};
    for(int n : cross_n) {
       std::string nstr = std::to_string(n);
       add("cross_r1.0_"+nstr, "cross polytope, dim-"+nstr+", oneNorm(x) <= 1", [n]() {
           Solved_Body* sb = generate_cross_polytope(n);
           sb->is_normalized = false;
           return sb;
       });
       add("cross_rn_"+nstr, "cross polytope, dim-"+nstr+", oneNorm(x) <= n [normalized]", [n]() {
           Solved_Body* s = generate_cross_polytope(n);
           Solved_Body* sb = s->scale(1.0/n);
           delete s;
           sb->is_normalized = true;
           return sb;
       });
       add("cross_rot_rn_"+nstr, "cross polytope, dim-"+nstr+", randomly rotated, oneNorm(x) <= n [normalized]", [n]() {
           Solved_Body* s1 = generate_cross_polytope(n);
           Solved_Body* s2 = s1->scale(1.0/n);
           Solved_Body* sb = s2->rotate();
           delete s1;
           delete s2;
           sb->is_normalized = true;
           return sb;
       });
    }

    // simplex_polytope
    std::vector<int> simplex_n = {2,3,10,20,40};
    for(int n : simplex_n) {
       std::string nstr = std::to_string(n);
       add("simplex_"+nstr, "simplex polytope, dim-"+nstr+", x_i>=0, oneNorm(x)<=2", [n]() {
           Solved_Body* sb = generate_simplex(n);
           sb->is_normalized = false;
           return sb;
       });
       add("simplex_preprocessed_"+nstr, "simplex polytope, dim-"+nstr+", x_i>=0, oneNorm(x)<=2 [preprocessed]", [n]() {
           Solved_Body* s = generate_simplex(n);
           Solved_Body* sb = s->preprocess();
           delete s;
           sb->is_normalized = true;
           return sb;
       });
    } 
    

    // ball - ellipsoids:
    std::vector<int> ball_n = {3,10,20,40,60,100,150,200};
    for(int n : ball_n) {
       std::string nstr = std::to_string(n);
       add("ball_r1.0_"+nstr, "ball ellipsoid, centered, dim-"+nstr+", radius 1 [normalized]", [n]() {
           Solved_Body* sb = generate_centered_ball(n,1.0);
           sb->is_normalized = true;
           return sb;
       });
       add("ball_rn_"+nstr, "ball ellipsoid, centered ,dim-"+nstr+", radius n [normalized]", [n]() {
           Solved_Body* s = generate_centered_ball(n,1.0);
           Solved_Body* sb = s->scale(1.0/n);
           delete s;
           sb->is_normalized = true;
           return sb;
       });
    }

    // random - ellipsoids:
    std::vector<int> ellipsoid_n = {3,10,20,40,60,100,150,200};
    for(int n : ball_n) {
       std::string nstr = std::to_string(n);
       add("ellipsoid_"+nstr, "random ellipsoid, centered, dim-"+nstr+", radius randomized", [n]() {
           Solved_Body* sb = generate_randomized_ellipsoid(n);
           sb->is_normalized = true;
           return sb;
       });
    }

    // half-ball / ellipsoids
    std::vector<int> half_n = {2,3,4,5,10,20,40,60,100};
    for(int n : half_n) {
       std::string nstr = std::to_string(n);
       add("half_preprocessed_"+nstr, "half-ball (ball+cube), dim-"+nstr+" [normalized]", [n]() {
           Solved_Body* ball = generate_centered_ball(n,1.0);
           FT lb[n];
           FT ub[n];
	   for(int i=0;i<n;i++) {lb[i] = -1; ub[i]=1;}
	   lb[0]=0; // truncate
	   ub[0]=1;
	   Solved_Body* c0 = generate_hyperrectangle(n, lb,ub);
           
	   Solved_Body* h0 = c0->join(ball);
	   h0->volume = ball->volume / 2.0;
	   
	   Solved_Body* sb = h0->preprocess();
           
	   delete ball;
	   delete c0;
	   delete h0;
	   return sb;
       });
       add("half_"+nstr, "half-ball (ball+cube), dim-"+nstr, [n]() {
           Solved_Body* ball = generate_centered_ball(n,1.0);
           FT lb[n];
           FT ub[n];
	   for(int i=0;i<n;i++) {lb[i] = -1; ub[i]=1;}
	   lb[0]=0; // truncate
	   ub[0]=1;
	   Solved_Body* c0 = generate_hyperrectangle(n, lb,ub);
           
	   Solved_Body* h0 = c0->join(ball);
	   h0->volume = ball->volume / 2.0;
	   
	   delete ball;
	   delete c0;
	   return h0;
       });
    }

    // 2-sphere
    std::vector<int> twosphere_n = {3,10,20,40,60,100,150,200,250,300};
    for(int n : twosphere_n) {
       std::string nstr = std::to_string(n);
       add("2sphere_preprocessed_"+nstr, "2 spheres, dim-"+nstr+" [normalized]", [n]() {
           Solved_Body* c0 = generate_centered_ball(n,1.0);
           Solved_Body* c2 = generate_centered_ball(n,1.0);
	   
           FT* a = (FT*)(aligned_alloc(32, n*sizeof(FT))); // align this to 32
	   for(int i=0;i<n;i++) {a[i]=(i==0);}
	   Solved_Body* c1 = c2->translate(a);
	   
	   Solved_Body* h0 = c0->join(c1);
	   h0->volume = 0;

	   // can calculate the volume as two n-ball segments
	   //
	   // S(n+1) = integrate Vn * sqrt(1-x^2)^n dx from 0.5 to 1
	   //        = Vn * 2 * integrate (1-t)^(n/2)*t^(-1/2) 1/2 dt from 0.25 to 1
	   //        = Vn * integrate (1-t)^(n/2)*t^(-1/2) dt from 0.25 to 1
	   // 
	   // calculate integral online
	   //
	   // links:
	   // https://math.stackexchange.com/questions/15656/volumes-of-n-balls-what-is-so-special-about-n-5
	   // https://en.wikipedia.org/wiki/Volume_of_an_n-ball
	   // https://www.integral-calculator.com/

	   std::map<int,FT> integrals = {
	      {2,  0.6141848493043784},
	      {3,  0.4166666666666667},
	      {4,  0.2982588737687016},
	      {5,  0.2208333333333333},
	      {10, 0.06329139154289489}, // up to here accurate
	      {20, 0.008368361860095765}, // seems to get worse from here...
	      {40, 0.0002521299865559373},
	      {60, 9.713735947843199 * 1e-6},
	      {100,1.890073619906658 * 1e-8},
	   };
	   auto it = integrals.find(n);
	   if(it != integrals.end()) {
	      FT ball = Ball_volume(n-1,1.0);
	      h0->volume = ball * it->second;
	   }
	   
	   Solved_Body* sb = h0->preprocess();
           
	   delete c0;
	   delete c1;
	   delete c2;
	   delete h0;
	   free(a);
	   return sb;
       });
    }

    // 2n-sphere
    std::vector<int> twonsphere_n = {2,3,4,5};
    for(int n : twonsphere_n) {
       std::string nstr = std::to_string(n);
       add("2nsphere_preprocessed_"+nstr, "2*n spheres, dim-"+nstr+" [normalized]", [n]() {
           Solved_Body* c0 = NULL;//generate_centered_ball(n,1.0);
           FT* a = (FT*)(aligned_alloc(32, n*sizeof(FT))); // align this to 32
	   
	   for(int j=0;j<n;j++) {
	      Solved_Body* c1 = generate_centered_ball(n,1.0);
	      Solved_Body* c2 = generate_centered_ball(n,1.0);
	      
	      for(int i=0;i<n;i++) {a[i]=0;}
	      
	      a[j] = 0.6;
	      Solved_Body* c1t = c1->translate(a);
	      a[j] = -0.6;
	      Solved_Body* c2t = c2->translate(a);
	      
	      Solved_Body* c3 = c1t->join(c2t);
	      Solved_Body* tmp = c0;

	      if(c0==NULL) {
	         c0 = c3;
	      } else {
	         c0 = c0->join(c3);
	         delete c3;
	         delete tmp;
	      }
	      delete c1;
	      delete c1t;
	      delete c2;
	      delete c2t;
	   }

	   Solved_Body* s0 = c0->shear();

	   Solved_Body* sb = s0->preprocess();
           sb->volume = 0;

	   delete s0;
	   delete c0;
	   free(a);
	   return sb;
       });
    }

    // 2-box
    std::vector<int> twobox_n = {3,10,20,40,60,100};
    for(int n : twobox_n) {
       std::string nstr = std::to_string(n);
       add("2box_"+nstr, "2 boxes, dim-"+nstr+"", [n]() {
           FT lb[n];
           FT ub[n];
	   for(int i=0;i<n;i++) {lb[i] = -1; ub[i]=1;}
	   lb[0]=0; // truncate
	   ub[0]=1;
	   Solved_Body* c0 = generate_hyperrectangle(n, lb,ub);
	   for(int i=0;i<n;i++) {lb[i] = -1; ub[i]=1;}
	   lb[1]=0; // truncate
	   ub[1]=1;
	   Solved_Body* c1 = generate_hyperrectangle(n, lb,ub);
           
	   Solved_Body* c2 = c0->join(c1);
	   c2->volume = c1->volume / 2.0;
	   Solved_Body* sb = c2->preprocess();
           
	   delete c0;
	   delete c1;
	   delete c2;
	   return sb;
       });
    }


    // Polyvest Polytopes:
    std::vector<std::string> polyvest_list = {
        "cc_8_10",
        "cc_8_11",
        "cross_13",
        "cross_7",
        "cross_9",
        "cube_10",
        "cube_10_2",
        "cube_14",
        "cube_14_2",
        "cube_15",
        "cube_2",
        "cube_20",
        "cube_25",
        "cube_30",
        "cube_35",
        "cube_40",
        "cube_5",
        "cube_80",
        "ex_1",
        "ex_2",
        "fm_6",
        "rect_3",
        "rh_1",
        "rh_2",
        "rh_20_40",
        "rh_3",
        "rh_30_60",
        "rh_4",
        "rh_40_80",
        "simplex_10",
        "simplex_14",
        "simplex_15",
        "simplex_20",
        "gagandeep"
    };
    for(auto f : polyvest_list) {
        std::string pname = "polyvest_"+f;
	add(pname, "One of polyvest polytopes, read from file.",[f]() {
	    if(const char* env_p = std::getenv("POLYVEST_PATH")) {
	        std::string path = env_p;
	        std::string fpath = path+f;
                //
                return generate_read_polyvest_polytope(fpath);
            } else {
	        std::cout << "ERROR: POLYVEST_PATH not set!\n" << std::endl;
		std::cout << "try: export POLYVEST_PATH='../polyvest/examples/'\n";
		std::exit(0);
	    }
	});
    }

    // vinci polytopes
    std::vector<std::string> vinci_list = {
                                           "cc_8_10.ine",
                                           "cc_8_11.ine",
                                           "cc_8_5.ine",
                                           "cc_8_6.ine",
                                           "cc_8_7.ine",
                                           "cc_8_8.ine",
                                           "cc_8_9.ine",
                                           "ccp_5.ine",
                                           "ccp_6.ine",
                                           /*"ccp_7.ine", doesn't have volume, took too long to compute...*/
                                           "Fm_4.ine",
                                           "Fm_5.ine",
                                           "Fm_6.ine",
                                           "rh_10_20.ine",
                                           "rh_10_25.ine",
                                           "rh_10_30.ine",
                                           "rh_8_20.ine",
                                           "rh_8_25.ine",
                                           "rh_8_30.ine",
                                           "rv_10_12.ine",
                                           "rv_10_13.ine",
                                           "rv_10_14.ine",
                                           "rv_8_10.ine",
                                           "rv_8_11.ine",
                                           "rv_8_12.ine",
                                           "rv_8_13.ine",
                                           "rv_8_14.ine",
                                           "rv_8_20.ine",
                                           "rv_8_30.ine",
                                           "rand_3_30_1000.ine",
                                           "rand_3_30_500.ine",
                                           "rand_3_30_100.ine",
                                           "rand_5_30_20.ine",
                                           "rand_7_30_20.ine",
                                           "rand_8_30_20.ine",

    };
    for(auto f : vinci_list) {
        std::string pname = "vinci_"+f;
	add(pname, "One of vinci polytopes, read from file.",[f]() {
	    if(const char* env_p = std::getenv("POLYVEST_PATH")) {
	        std::string path = env_p;
	        std::string fpath = path+f;
                //
                return generate_read_vinci_polytope(fpath);
            } else {
	        std::cout << "ERROR: POLYVEST_PATH not set!\n" << std::endl;
		std::cout << "try: export POLYVEST_PATH='../polyvest/examples/'\n";
		std::exit(0);
	    }
	});
    }

    for(int n=3;n<=20;n++) {
	std::string nstr = std::to_string(n);
        std::string pname = "birk_"+nstr;
	add(pname, "One of birk polytopes, from volumeEsti, read from file.",[nstr]() {
	    if(const char* env_p = std::getenv("POLYVEST_PATH")) {
	        std::string path = env_p;
	        std::string fpath = path+"/birk/birk"+nstr+".ine";
                //
                return generate_read_vinci_polytope(fpath);
            } else {
	        std::cout << "ERROR: POLYVEST_PATH not set!\n" << std::endl;
		std::cout << "try: export POLYVEST_PATH='../polyvest/examples/'\n";
		std::exit(0);
	    }
	});
    }


    // kvariable
    std::vector<int> kvar_n = {4,5,10,20,30,40,50,60,100,150,200};
    for(int n : kvar_n) {
       std::string nstr = std::to_string(n);
       add("2var_TSP_"+nstr, "2-variable-polytope, translated, axisScaled, preprocessed, 10n constraints, "+nstr+"-dim [normalized]", [n]() {
           Solved_Body* b1 = generate_kvariable_polytope(n,2,1.0,10*n);//k=2, r=1.0
           
	   std::cout << "\n# Generator: b1:\n";
	   //b1->print();

	   FT* a = (FT*)(aligned_alloc(32, n*sizeof(FT))); // align this to 32
	   
	   for(int i=0;i<n;i++) {a[i]=prng_get_random_double_in_range(-10,10);}
	   Solved_Body* b2 = b1->translate(a);
	   
	   std::cout << "\n# Generator: b2 - translated:\n";
	   //b2->print();

	   for(int i=0;i<n;i++) {a[i]=prng_get_random_double_in_range(0.1,10);}
	   Solved_Body* b3 = b2->scaleAxis(a);
	   
	   std::cout << "\n# Generator: b3 - scaleAxis:\n";
	   //b3->print();
	   
	   Solved_Body* sb = b3->preprocess();
           
	   delete b1;
	   delete b2;
	   delete b3;
	   free(a);
           return sb;
       });
       add("xvar_TSR_"+nstr, "2-variable-polytope, translated, axisScaled, rotated, 10n constraints, "+nstr+"-dim [for preprocessing test]", [n]() {
           Solved_Body* b1 = generate_kvariable_polytope(n,2,1.0,10*n);//k=2, r=1.0
           
	   std::cout << "\n# Generator: b1:\n";
	   //b1->print();

	   FT* a = (FT*)(aligned_alloc(32, n*sizeof(FT))); // align this to 32
	   
	   for(int i=0;i<n;i++) {a[i]=prng_get_random_double_in_range(-10,10);}
	   Solved_Body* b2 = b1->translate(a);
	   
	   std::cout << "\n# Generator: b2 - translated:\n";
	   //b2->print();

	   for(int i=0;i<n;i++) {a[i]=std::pow(2.0,prng_get_random_double_in_range(-5,5));}
	   Solved_Body* b3 = b2->scaleAxis(a);
	   
	   std::cout << "\n# Generator: b3 - scaleAxis:\n";
	   //b3->print();
	   
	   Solved_Body* sb = b3->rotate();
           
	   delete b1;
	   delete b2;
	   delete b3;
	   free(a);
           return sb;
       });
 
       std::vector<std::string> TSP_precomputed = {"2var_TSP_100.ine",
                                                   "2var_TSP_200.ine",
                                                   "2var_TSP_40.ine",
                                                   "2var_TSP_10.ine",
                                                   "2var_TSP_20.ine",
                                                   "2var_TSP_50.ine",
                                                   "2var_TSP_150.ine",
                                                   "2var_TSP_30.ine",
                                                   "2var_TSP_60.ine",};
       
       for(auto f : TSP_precomputed) {
           std::string pname = "vinci_"+f;
           add(pname, "precomputed 2var_TSP, read from file.",[f]() {
                                                                    if(const char* env_p = std::getenv("POLYVEST_PATH")) {
                                                                        std::string path = env_p;
                                                                        std::string fpath = path+f;
                                                                        //
                                                                        Solved_Body *sb = generate_read_vinci_polytope(fpath);
                                                                        sb->is_normalized = true;
                                                                        return sb;
                                                                    } else {
                                                                        std::cout << "ERROR: POLYVEST_PATH not set!\n" << std::endl;
                                                                        std::cout << "try: export POLYVEST_PATH='../polyvest/examples/'\n";
                                                                        std::exit(0);
                                                                    }
                                                                });
       }
       add("2var_"+nstr, "2-variable polytope, 10n constraints,"+nstr+"-dim [normalized]", [n]() {
           Solved_Body* sb = generate_kvariable_polytope(n,2,1.0,10*n);//k=2, r=1.0
           sb->is_normalized = true;
           return sb;
       });
       add("3var_"+nstr, "3-variable polytope, 10n constraints,"+nstr+"-dim [normalized]", [n]() {
           Solved_Body* sb = generate_kvariable_polytope(n,3,1.0,10*n);//k=3, r=1.0
           sb->is_normalized = true;
           return sb;
       });
       add("4var_"+nstr, "4-variable polytope, 10n constraints,"+nstr+"-dim [normalized]", [n]() {
           Solved_Body* sb = generate_kvariable_polytope(n,4,1.0,10*n);//k=4, r=1.0
           sb->is_normalized = true;
           return sb;
       });
    }

    // control density
    std::vector<int> dens = {1,2,3,4,5,6,7,8,9,10,11};
    for (auto f : dens){
        std::string nstr = std::to_string(f);
        int density = min((int) pow(1.55,f), 100);
        add("dens_"+nstr,
            "100-dim polytope with density " + std::to_string(density) + " [normalized]",
            [density]()
            {
                Solved_Body* sb = generate_kvariable_polytope(100,density,1.0,1000,false);//k=2, r=1.0
                sb->is_normalized = true;
                return sb;
            }
            );

        add("dens200_"+nstr,
            "200-dim polytope with density " + std::to_string(density) + " [normalized]",
            [density]()
            {
                Solved_Body* sb = generate_kvariable_polytope(200,2*density,1.0,2000,false);//k=2, r=1.0
                sb->is_normalized = true;
                return sb;
            }
            );
        
    }

    
}

Solved_Body*
Solved_Body::clone() {
    Solved_Body* sb = new Solved_Body(bcount,n);
    sb->volume = volume;
    sb->is_normalized = is_normalized;
    for(int b=0; b<bcount; b++) {
       sb->type[b] = type[b];
       sb->body[b] = type[b]->clone(body[b]);
    }
    return sb;
}

Solved_Body*
Solved_Body::join(const Solved_Body* other) {
    assert(n==other->n);
    Solved_Body* sb = new Solved_Body(bcount + other->bcount,n);
    sb->volume = 0;// volume unknown
    sb->is_normalized = false;
    for(int b=0; b<bcount; b++) {
       sb->type[b] = type[b];
       sb->body[b] = type[b]->clone(body[b]);
    }
    for(int b=0; b<other->bcount; b++) {
       sb->type[bcount+b] = other->type[b];
       sb->body[bcount+b] = other->type[b]->clone(other->body[b]);
    }
    return sb;
}

Solved_Body*
Solved_Body::transform(const Matrix* L, const FT det, const FT* a, const FT beta) {
    Solved_Body* sb = clone();
    for(int b=0; b<bcount; b++) {
	type[b]->transform(body[b],sb->body[b],L,a,beta);
    }
    sb->volume = volume / std::pow(beta,n) / det;
    return sb;
}

Solved_Body*
Solved_Body::scale(const FT beta) {
    Matrix* L = Matrix_new(n,n);
    FT* a = (FT*)(aligned_alloc(32, n*sizeof(FT))); // align this to 32
    for(int i=0;i<n;i++) {
        Matrix_set(L, i, i, 1.0);
        a[i] = 0;
    }
    Solved_Body* sb = transform(L,1.0,a,beta);
    free(a);
    Matrix_free(L);
    return sb;
}


Solved_Body*
Solved_Body::scaleAxis(const FT* diag) {
    Matrix* L = Matrix_new(n,n);
    FT* a = (FT*)(aligned_alloc(32, n*sizeof(FT))); // align this to 32
    for(int i=0;i<n;i++) {
        Matrix_set(L, i, i, diag[i]);
        a[i] = 0;
    }
    Solved_Body* sb = transform(L,1.0,a,1.0);
    free(a);
    Matrix_free(L);
    return sb;
}



Solved_Body*
Solved_Body::rotate() {
    Matrix* L = Matrix_new(n,n);
    FT* a = (FT*)(aligned_alloc(32, n*sizeof(FT))); // align this to 32
    for(int i=0;i<n;i++) {
        Matrix_set(L, i, i, 1.0);
        a[i] = 0;
    }
    for(int i=0;i<n;i++) {
        for(int j=0;j<n;j++) {
	    FT angle = prng_get_random_double_in_range(0,2*M_PI); 
            Matrix_rotate(L, i, j, angle);
        }
    }
    for(int b=0; b<bcount; b++) {
       assert(type[b]==&Polytope_T || type[b]==&PolytopeT_T);
       // it is not safe to transform Ellipsoid with non-L matrix
    }
    Solved_Body* sb = transform(L,1.0,a,1.0);
    free(a);
    Matrix_free(L);
    return sb;
}

Solved_Body*
Solved_Body::shear() {
    Matrix* L = Matrix_new(n,n);
    FT* a = (FT*)(aligned_alloc(32, n*sizeof(FT))); // align this to 32
    for(int i=0;i<n;i++) {
	for(int j=0;j<i;j++) {
	   FT s = prng_get_random_double_in_range(-0.6,0.6); 
           Matrix_set(L, i, j, s);
	}
        Matrix_set(L, i, i, 1.0);
        a[i] = 0;
    }
    Solved_Body* sb = transform(L,1.0,a,1.0);
    free(a);
    Matrix_free(L);
    return sb;
}



Solved_Body*
Solved_Body::translate(const FT* a) {
    Matrix* L = Matrix_new(n,n);
    for(int i=0;i<n;i++) {
        Matrix_set(L, i, i, 1.0);
    }
    Solved_Body* sb = transform(L,1.0,a,1.0);
    Matrix_free(L);
    return sb;
}



Solved_Body*
Solved_Body::preprocess() {
    Solved_Body* sb = clone();
    ArbitraryExpNum proc_det = ArbitraryExpNum_new(1);
    preprocess_generic(n, bcount, (const void**)body, sb->body, (const Body_T**)type, &proc_det);
    sb->is_normalized = true;
    sb->volume = volume / proc_det.num;
    return sb;
}

void
Solved_Body::polytopeTranspose() {
    for(int b=0; b<bcount; b++) {
       if(type[b]==&Polytope_T) {
          Polytope* p = (Polytope*)body[b];
	  PolytopeT* pt = Polytope_to_PolytopeT(p);
	  Polytope_free(p);
	  body[b] = pt;
	  type[b] = &PolytopeT_T;
       }
    }
}

void
Solved_Body::polytopeCSC(){
    for (int b = 0; b < bcount; b++){
        if (type[b] == &Polytope_T){
            Polytope *p = (Polytope *) body[b];
            PolytopeCSC *pcsc = Polytope_to_PolytopeCSC(p);
            Polytope_free(p);
            body[b] = pcsc;
            type[b] = &PolytopeCSC_T;
        }
    }
}

void
Solved_Body::polytopeJIT() {
    for(int b=0; b<bcount; b++) {
       if(type[b]==&Polytope_T) {
          Polytope* p = (Polytope*)body[b];
	  PolytopeJIT* pjit = Polytope_to_PolytopeJIT(p);
	  Polytope_free(p);
	  body[b] = pjit;
	  type[b] = &PolytopeJIT_T;
       }
    }
}

void
Solved_Body::optimize() {
    assert(bcount == 1 && "only one body allowed for Polytope optimize");
    assert(type[0] == &Polytope_T && "body must be Polytope");
    Polytope* p = (Polytope*)body[0];
    Polytope* q = optimize_polytope(p);
    Polytope_free(p);
    body[0] = q;
}




Solved_Body_Generator* solved_body_generator_ = NULL;
Solved_Body_Generator* solved_body_generator() {
    if(!solved_body_generator_) {solved_body_generator_ = new Solved_Body_Generator();}
    return solved_body_generator_;
}

Solved_Body* generate_hyperrectangle(int dims, FT *lower_bounds, FT *upper_bounds) {
    int num_constraints = 2*dims;

    Polytope *hyperrectangle = Polytope_new(dims, num_constraints);

    // Lower bounds
    for (int i = 0; i < dims; i++) {
        // Ax >= b means that -Ax <= -b
        Polytope_set_a(hyperrectangle, i, i, -1);
        Polytope_set_b(hyperrectangle, i, -lower_bounds[i]);
    }

    // Upper bounds
    for (int i = 0; i < dims; i++) {
        Polytope_set_a(hyperrectangle, i + dims, i, 1);
        Polytope_set_b(hyperrectangle, i + dims, upper_bounds[i]);
    }

    FT volume = 1.0;
    for (int i = 0; i < dims; i++) {
        volume *= upper_bounds[i] - lower_bounds[i];
    }

    Solved_Body* result = new Solved_Body(1,dims);
    result->body[0] = hyperrectangle;
    result->type[0] = &Polytope_T;
    result->volume = volume;
    return result;
}

Solved_Body* generate_centered_hypercube(int dims, FT r) {
    FT lower_bounds[dims];
    FT upper_bounds[dims];

    for (int i = 0; i < dims; i++) {
        lower_bounds[i] = -r;
        upper_bounds[i] = +r;
    }

    return generate_hyperrectangle(dims, lower_bounds, upper_bounds);
}

Solved_Body* generate_cross_polytope(int dims) {

    // This cross-polytope is defined by ± x_1 ± x_2 ... ± x_n <= 1
    int num_constraints = (1 << dims);

    Polytope *cross_polytope = Polytope_new(dims, num_constraints);

    for (int i = 0; i < num_constraints; i++) {
        for (int j = 0; j < dims; j++) {
            // i is a number from 0 to 2^dims - 1
            // We check its j-th bit to decide whether we want x_j or -x_j
            FT plus_minus_1 = 1.0 - 2 * ((i & (1 << j)) >> j);
            Polytope_set_a(cross_polytope, i, j, plus_minus_1);
        }
    }

    for (int i = 0; i < num_constraints; i++) {
        Polytope_set_b(cross_polytope, i, 1);
    }

    // Volume of a cross polytope is 2^n / n!
    long n_factorial = 1;
    for (long i = 1; i <= dims; i++) {
        n_factorial *= i;
    }
    FT volume = ((FT) num_constraints) / n_factorial;

    int bcount = 1;
    Solved_Body *result = new Solved_Body(bcount, dims);
    result->body[0] = cross_polytope;
    result->type[0] = &Polytope_T;
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

    Polytope *simplex = Polytope_new(dims, num_constraints);

    // x_i >= 0
    for (int i = 0; i < dims; i++) {
        Polytope_set_a(simplex, i, i, -1);
        Polytope_set_b(simplex, i, 0);
    }

    // x_1 + ... x_n <= 2
    for (int i = 0; i < dims; i++) {
        Polytope_set_a(simplex, dims, i, 1);
    }
    Polytope_set_b(simplex, dims, 2);

    // Again, volume = 2^n / n!, just like the cross polytope!
    long n_factorial = 1;
    for (long i = 1; i <= dims; i++) {
        n_factorial *= i;
    }
    FT volume = std::pow(2,dims) / n_factorial;

    int bcount = 1;
    Solved_Body *result = new Solved_Body(bcount, dims);
    result->body[0] = simplex;
    result->type[0] = &Polytope_T;
    result->volume = volume;
    return result;

}

struct MyElement {
    int index;
    double r;

    bool operator > (const MyElement& other) const {
       return r > other.r;
    }
};

void choosek(const int n, const int k, std::vector<int> &choice) {
    std::vector<MyElement> e;
    e.reserve(n);
    for(int i=0;i<n;i++) {
        e.push_back({i,prng_get_random_double_in_range(0,1)});
    }
    std::sort(e.begin(),e.end(),greater<MyElement>());
    choice.resize(k,0);
    for(int i=0;i<k;i++) {choice[i] = e[i].index;}
}

Solved_Body* generate_kvariable_polytope(const int dims, const int k, const FT r, const int num_constraints, const bool boundingBox) {
    //assert(k>=2);
    assert(num_constraints >= 2*dims);
    Polytope *p = Polytope_new(dims, num_constraints);
    
    // add cube at the end to make sure polytope is bounded.
    const int rand_constr = num_constraints - boundingBox*2*dims;

    int j = 0;
    std::vector<int> choice;
    FT* d = (FT*)(aligned_alloc(32, dims*sizeof(FT))); // align this to 32
    for(int c=0;c<rand_constr;c++) {
        choosek(dims,k,choice);
	for(int i=0;i<dims;i++) {d[i]=0;}
	for(int i=0;i<k;i++) {d[choice[i]] = prng_get_random_double_normal();}
	FT d2 = squaredNorm(d,dims);
	FT dd = sqrt(d2);
	for(int i=0;i<k;i++) {Polytope_set_a(p, j, choice[i], d[choice[i]]/dd);}
	Polytope_set_b(p, j, r);
	j++;
    }
    free(d);
    
    if(boundingBox) {
       // add the cube:
       for(int i=0;i<dims;i++) {
           Polytope_set_a(p, j, i, -1);
           Polytope_set_b(p, j, r);
           j++;
           Polytope_set_a(p, j, i, 1);
           Polytope_set_b(p, j, r);
           j++;
       }
    }
    assert(j==num_constraints);
    
    FT volume = std::pow(2,dims);

    int bcount = 1;
    Solved_Body *result = new Solved_Body(bcount, dims);
    result->body[0] = p;
    result->type[0] = &Polytope_T;
    result->volume = volume;
    return result;
}

Solved_Body* generate_ellipsoid(int dims, FT *lower_bounds, FT *upper_bounds) {

    Ellipsoid *ellipsoid = Ellipsoid_new_with_T(dims);

    FT volume = Ball_volume(dims,1.0);// unit ball volume

    for (int i = 0; i < dims; i++) {
        FT radius_i = (upper_bounds[i] - lower_bounds[i]) / 2;
        FT midpoint = lower_bounds[i] + radius_i;

        ellipsoid->a[i] = midpoint;
        
        // Instead of x_1^2 + ... + x_n^2 <= 1
        // we want that (x_1/r_1)^2 + ... + (x_n/r_n)^2 <= 1
        FT lambda_i = 1 / (radius_i * radius_i);
        Ellipsoid_set_a(ellipsoid, i, i, lambda_i);
        Ellipsoid_set_Ta(ellipsoid, i, i, lambda_i);

        // Over all i's we calculate det(T) here, where T = A^{-1}
        volume *= radius_i * radius_i;
    }

    int bcount = 1;
    Solved_Body *result = new Solved_Body(bcount, dims);
    result->body[0] = ellipsoid;
    result->type[0] = &Ellipsoid_T;
    result->volume = volume;
    return result;

}

Solved_Body* generate_centered_ball(int dims, FT r) {

    FT lower_bounds[dims];
    FT upper_bounds[dims];

    for (int i = 0; i < dims; i++) {
        lower_bounds[i] = -r;
        upper_bounds[i] = +r;
    }

    return generate_ellipsoid(dims, lower_bounds, upper_bounds);

}

Solved_Body* generate_randomized_ellipsoid(int dims) {

    Ellipsoid *ellipsoid = Ellipsoid_new_with_T(dims);

    FT volume = Ball_volume(dims,1.0);// unit ball volume

    for (int i = 0; i < dims; i++) {
        ellipsoid->a[i] = 0.0;
        FT* Ai = Ellipsoid_get_Ai(ellipsoid,i);
        Ai[i] = prng_get_random_double_in_range(1.0,10.0);
        Ellipsoid_set_Ta(ellipsoid, i, i, Ai[i]);
    }

    int bcount = 1;
    Solved_Body *result = new Solved_Body(bcount, dims);
    result->body[0] = ellipsoid;
    result->type[0] = &Ellipsoid_T;
    result->volume = volume;
    return result;

}

Solved_Body* generate_read_polyvest_polytope(const std::string &fileName) {
    Polytope *P;
    int err = read_polyvest_p(fileName, &P);
    assert(!err &&
           "couldn't read example polytope");

    Solved_Body *result = new Solved_Body(1, P->n);
    result->body[0] = P;
    result->type[0] = &Polytope_T;
    result->volume = 0; // unknown
    
    return result;
}

Solved_Body *generate_read_vinci_polytope(const std::string &filename){
    Polytope *P;
    FT vol;
    int err = read_vinci(filename, &P, &vol);
    assert(!err && "couldn't read vinci polytope");
    Solved_Body *res = new Solved_Body(1, P->n);
    res->body[0] = P;
    res->type[0] = &Polytope_T;
    res->volume = vol;

    return res;
}


//
//// This function is out of order, because read_polyvest_p() doesn't exist anymore
//// We need to reimplement it anyway to return Polytope instead of Polytope
///*
//struct Solved_Body generate_solved_polyvest_polytope(int index) {
//
//    Polytope *polytope;
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

