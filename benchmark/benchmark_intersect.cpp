#include "benchmark.hpp"
#include <immintrin.h>
#include "../src/volume/volume_helper.hpp"


class Benchmark_intersect : public Benchmark_base {
    public:
        Benchmark_intersect(std::string name, int reps, bool convergence, int warmup_reps, const std::string &generator, const int polytopeType, const bool polytopeOptimize, const double time_ci_alpha_, const double results_ci_alpha_, const bool printBody)
		: Benchmark_base(name, reps, convergence, warmup_reps, time_ci_alpha_, results_ci_alpha_), generator(generator), polytopeType(polytopeType), polytopeOptimize(polytopeOptimize), printBody(printBody) {}

    protected:
        void initialize () {
            std::cout << "initializing intersect data..." << std::endl;
        prng_init();
	   
            solved_body = solved_body_generator()->get(generator,false);
	    if(polytopeOptimize) {
	       solved_body->optimize();
	    }

	    switch(polytopeType) {
            case 0: // row major
                break;
            case 1: // column major
                solved_body->polytopeTranspose();
                break;
            case 2: // CSC format
                solved_body->polytopeCSC();
                break;
            case 3: // JIT format
                solved_body->polytopeJIT();
                break;
	    }
	    if(printBody) {solved_body->print();}
	    assert(solved_body->is_normalized);
	    assert(solved_body->bcount == 1);
	    
	    int n = solved_body->n; // allocate more memory for 4/8 sets!
	    x = (FT*)(aligned_alloc(32, 8*n*sizeof(FT))); // align this to 32
	    d = (FT*)(aligned_alloc(32, 8*n*sizeof(FT))); // align this to 32
	    
	    int cache_size = solved_body->type[0]->cacheAlloc(solved_body->body[0]);
            cache = aligned_alloc(32, 8*cache_size); // align this to 32

	    reset();
        }
        void reset () {
	    int n = solved_body->n;
            for(int i=0; i<n;i++) {
	        x[i] = prng_get_random_double_0_1()*0.1;
	        d[i] = prng_get_random_double_normal();
	    }
	    solved_body->type[0]->cacheReset(solved_body->body[0], x, cache);
	    dd = prng_get_random_int_in_range(0,n-1);
	}
        double run () {
	    int n = solved_body->n;
	    FT t0, t1;
	    solved_body->type[0]->intersect(solved_body->body[0], x, d, &t0, &t1);
            
	    // step now
	    FT t = prng_get_random_double_in_range(t0,t1);
            for(int j=0;j<n;j++) {x[j] += d[j]*t;}
            return 0;
	}
	void finalize() {
	    int n = solved_body->n;
	    pc_stack().reset();
            {
		{
		   PC_Frame<intersect_cost_f> frame((void*)solved_body->type[0]->intersect);
                   frame.costf()(solved_body->body[0]);
                }
                //pc_stack().log(0,0, "random double - TODO");
                {// frame for random double_in_range
                    PC_Frame<random_double_in_range_cost_f> frame((void*) prng_get_random_double_in_range);
                    frame.costf()(NULL);
                }        
                pc_stack().log(2*n, 3*n*sizeof(FT)," x += d*t");
	    }
            pc_stack().print();
	    pc_flops = pc_stack().flops();
	    pc_bytes = pc_stack().bytes();
	}
    protected:
	const std::string generator;
	Solved_Body* solved_body;
	FT* x;
	FT* d;
	int dd;
	void* cache;
	int polytopeType = 0;
	const bool polytopeOptimize;
	const bool printBody;
};

class Benchmark_intersectCoord : public Benchmark_intersect {
    public:
        Benchmark_intersectCoord(std::string name, int reps, bool convergence, int warmup_reps, const std::string &generator, const int polytopeType, const bool polytopeOptimize, const double time_ci_alpha_, const double results_ci_alpha_, const bool printBody, const bool intersectCoord_intersect, const bool intersectCoord_update)
		: Benchmark_intersect(name, reps, convergence, warmup_reps, generator, polytopeType, polytopeOptimize, time_ci_alpha_, results_ci_alpha_, printBody), intersectCoord_intersect(intersectCoord_intersect), intersectCoord_update(intersectCoord_update) {}
    	void finalize() {
	    pc_stack().reset();
            {
		if(intersectCoord_intersect){
		    PC_Frame<intersect_cost_f> frame((void*)solved_body->type[0]->intersectCoord);
                    frame.costf()(solved_body->body[0]);
		}

                if(intersectCoord_update) {
		   //pc_stack().log(0, 0, "random double - TODO");
                {// frame for random double_in_range
                    PC_Frame<random_double_in_range_cost_f> frame((void*) prng_get_random_double_in_range);
                    frame.costf()(NULL);
                } 
	           // Reading and writing x[dd] with one add in between
                   pc_stack().log(1, 2, "x[dd] += t;");
                   
                   // body intersectCoord
		   {
	               PC_Frame<cacheUpdateCoord_cost_f> frame((void*) solved_body->type[0]->cacheUpdateCoord);
                       frame.costf()(solved_body->body[0]);
		   }
		}
	    }
            pc_stack().print();
	    pc_flops = pc_stack().flops();
	    pc_bytes = pc_stack().bytes();
	}
    protected:
        double run () {
	    FT t0=-1, t1=1;
	    if(intersectCoord_intersect) {
	       solved_body->type[0]->intersectCoord(solved_body->body[0], x, dd, &t0, &t1, cache);
	    }
	    if(intersectCoord_update) {
	       // step now
	       FT t = prng_get_random_double_in_range(t0,t1);
               x[dd] += t;
               solved_body->type[0]->cacheUpdateCoord(solved_body->body[0], dd, t, cache);
	    }
	    return 0;
	}
	const bool intersectCoord_intersect;
	const bool intersectCoord_update;
};


class Benchmark_intersectCoord4 : public Benchmark_intersectCoord {
    public:
        Benchmark_intersectCoord4(std::string name, int reps, bool convergence, int warmup_reps, const std::string &generator, const int polytopeType, const bool polytopeOptimize, const double time_ci_alpha_, const double results_ci_alpha_, const bool printBody, const bool intersectCoord_intersect, const bool intersectCoord_update)
		: Benchmark_intersectCoord(name, reps, convergence, warmup_reps, generator, polytopeType, polytopeOptimize, time_ci_alpha_, results_ci_alpha_, printBody,intersectCoord_intersect,intersectCoord_update) {}
	void reset () {
	    int n = solved_body->n;
            for(int i=0; i<4*n;i++) {
	        x[i] = prng_get_random_double_0_1()*0.1;
	        d[i] = prng_get_random_double_normal();
	    }
	    solved_body->type[0]->cacheReset4(solved_body->body[0], x, cache);
	    dd = prng_get_random_int_in_range(0,n-1);
	}
    
    	void finalize() {
	    pc_stack().reset();
            {
		if(intersectCoord_intersect){
		    PC_Frame<intersect_cost_f> frame((void*)solved_body->type[0]->intersectCoord4);
                    frame.costf()(solved_body->body[0]);
		}

                if(intersectCoord_update) {
                   {// frame for random double_in_range
                       PC_Frame<rand256d_cost_f_t> frame((void*) rand256d_f);
                       frame.costf()();
                   } 
	           // Reading and writing x[dd] with one add in between
                   pc_stack().log(4, 8*sizeof(FT), "x[dd] += t;");
                   
                   // body intersectCoord
		   {
	               PC_Frame<cacheUpdateCoord_cost_f> frame((void*) solved_body->type[0]->cacheUpdateCoord4);
                       frame.costf()(solved_body->body[0]);
		   }
		}
	    }
            pc_stack().print();
	    pc_flops = pc_stack().flops();
	    pc_bytes = pc_stack().bytes();
	}
    protected:
        double run () {
	    //FT t0=-1, t1=1;
	    __m256d lo = _mm256_set1_pd(0);
	    __m256d hi = _mm256_set1_pd(0);
	    if(intersectCoord_intersect) {
	       FTpair4 tp = solved_body->type[0]->intersectCoord4(solved_body->body[0], x, dd, cache);
	       lo = tp.low0;
	       hi = tp.hi0;
	    }
	    if(intersectCoord_update) {
	       // step now
	       __m256d t = rand256d_f();
	       t = _mm256_fmadd_pd(_mm256_sub_pd(hi,lo), t, lo);
	       
               //x[dd] += t;
	       __m256d xdd = _mm256_load_pd(x+dd*4);
	       __m256d xdd_t = _mm256_add_pd(xdd,t);
	       _mm256_store_pd(x+dd*4, xdd_t);
               solved_body->type[0]->cacheUpdateCoord4(solved_body->body[0], dd, t, cache);
	    }
	    return 0;
	}
};

class Benchmark_intersectCoord8 : public Benchmark_intersectCoord {
    public:
        Benchmark_intersectCoord8(std::string name, int reps, bool convergence, int warmup_reps, const std::string &generator, const int polytopeType, const bool polytopeOptimize, const double time_ci_alpha_, const double results_ci_alpha_, const bool printBody, const bool intersectCoord_intersect, const bool intersectCoord_update)
		: Benchmark_intersectCoord(name, reps, convergence, warmup_reps, generator, polytopeType, polytopeOptimize, time_ci_alpha_, results_ci_alpha_, printBody,intersectCoord_intersect,intersectCoord_update) {}
	void reset () {
	    int n = solved_body->n;
            for(int i=0; i<8*n;i++) {
	        x[i] = prng_get_random_double_0_1()*0.1;
	        d[i] = prng_get_random_double_normal();
	    }
	    solved_body->type[0]->cacheReset8(solved_body->body[0], x, cache);
	    dd = prng_get_random_int_in_range(0,n-1);
	}
    
    	void finalize() {
	    pc_stack().reset();
            {
		if(intersectCoord_intersect){
		    PC_Frame<intersect_cost_f> frame((void*)solved_body->type[0]->intersectCoord8);
                    frame.costf()(solved_body->body[0]);
		}

                if(intersectCoord_update) {
                   {// frame for random double_in_range
                       PC_Frame<rand256d_cost_f_t> frame((void*) rand256d_f,2);
                       frame.costf()();
                   } 
	           // Reading and writing x[dd] with one add in between
                   pc_stack().log(8, 16*sizeof(FT), "x[dd] += t;");
                   
                   // body intersectCoord
		   {
	               PC_Frame<cacheUpdateCoord_cost_f> frame((void*) solved_body->type[0]->cacheUpdateCoord8);
                       frame.costf()(solved_body->body[0]);
		   }
		}
	    }
            pc_stack().print();
	    pc_flops = pc_stack().flops();
	    pc_bytes = pc_stack().bytes();
	}
    protected:
        double run () {
	    //FT t0=-1, t1=1;
	    __m256d lo0 = _mm256_set1_pd(0);
	    __m256d lo1 = _mm256_set1_pd(0);
	    __m256d hi0 = _mm256_set1_pd(0);
	    __m256d hi1 = _mm256_set1_pd(0);
	    if(intersectCoord_intersect) {
	       FTpair8 tp = solved_body->type[0]->intersectCoord8(solved_body->body[0], x, dd, cache);
	       lo0 = tp.low0;
	       lo1 = tp.low1;
	       hi0 = tp.hi0;
	       hi1 = tp.hi1;
	    }
	    if(intersectCoord_update) {
	       // step now
	       __m256d t0 = rand256d_f();
	       t0 = _mm256_fmadd_pd(_mm256_sub_pd(hi0,lo0), t0, lo0);
	       __m256d t1 = rand256d_f();
	       t1 = _mm256_fmadd_pd(_mm256_sub_pd(hi1,lo1), t1, lo1);
	       //x[dd] += t;
	       __m256d xdd0 = _mm256_load_pd(x+dd*8);
	       __m256d xdd1 = _mm256_load_pd(x+dd*8+4);
	       __m256d xdd_t0 = _mm256_add_pd(xdd0,t0);
	       __m256d xdd_t1 = _mm256_add_pd(xdd1,t1);
	       _mm256_store_pd(x+dd*8,   xdd_t0);
	       _mm256_store_pd(x+dd*8+4, xdd_t1);
               solved_body->type[0]->cacheUpdateCoord8(solved_body->body[0], dd, {t0,t1}, cache);
	    }
	    return 0;
	}
};



int main(int argc, char *argv[]){
    CLI cli(argc,argv,"benchmark");
    CLIFunctionsVolume cliFun(cli);
    
    int r = 100;
    int warmup = 0;
    double time_ci_alpha;
    double results_ci_alpha;
    cliFun.claimOpt('b',"Benchmarking configuration");
    cliFun.add(new CLIF_OptionNumber<int>(&r,'b',"r","100", 1, 10000000));
    cliFun.add(new CLIF_OptionNumber<int>(&warmup,'b',"warmup","0", 0, 10000000));
    cliFun.add(new CLIF_OptionNumber<double>(&time_ci_alpha,'b',"time_ci_alpha","0.95", 0, 1));
    cliFun.add(new CLIF_OptionNumber<double>(&results_ci_alpha,'b',"results_ci_alpha","0.95", 0, 1));
    
    std::string generator = "cube";
    auto &gen_map = solved_body_generator()->gen_map();
    cliFun.add(new CLIF_Option<std::string>(&generator,'b',"generator","cube_r1.0_10", gen_map));
    
    
    bool intersectCoord = false;
    bool intersectCoord_intersect = true;
    bool intersectCoord_update = true;
    cliFun.add(new CLIF_TrippleOption<bool,bool,bool>(
			    &intersectCoord, &intersectCoord_intersect, &intersectCoord_update,
			    'b',"intersect","intersect", {
                                                     {"intersect",      {{false,{false,false}},  "random direction intersection"}},
                                                     {"intersectCoord", {{true ,{true ,true }},  "coordinate direction intersect + update"}},
                                                     {"intersectCoord_only", {{true ,{true ,false }},  "coordinate direction intersect only, no update"}},
                                                     {"cacheUpdateCoord", {{true ,{false ,true }},  "coordinate direction cache update, no intersect"}},
    }));
    
    int intersectSet = 1;
    cliFun.add(new CLIF_Option<int>(&intersectSet,'b',"intersectSet","1", {
                                                     {"1",{1,"single lane x, the conventional impl"}},
						     {"4", {4, "4-way x"}}, 
						     {"8", {8, "8-way x"}}, 
						     }));



    bool polytopeOptimize = false;
    cliFun.add(new CLIF_Option<bool>(&polytopeOptimize,'b',"polytopeOptimize","false", {
                                                     {"false",{false,"-"}},
						     {"true", {true, "Sort constraints to optimize access pattern"}} }));

    bool printBody = false;
    cliFun.add(new CLIF_Option<bool>(&printBody,'b',"printBody","false", {
                                                     {"false",{false,"-"}},
						     {"true", {true, "Print body before benchmark is run."}} }));



    int polytopeType = 0;
    cliFun.add(new CLIF_Option<int>(&polytopeType,'b',"polytopeType","0",
                                    {
                                     {"0",{0, "Polytope format / rows"}},
                                     {"1",{1, "PolytopeT format / columns"}},
                                     {"2",{2, "PolytopeCSC format"}},
                                     {"3",{3, "PolytopeJIT format"}},
                                    }));

    cliFun.preParse();
    if (!cli.parse()) {return -1;}
    cliFun.postParse();
    
    if(!intersectCoord) {
	assert(intersectSet == 1);
        Benchmark_intersect b("intersect", r, true, warmup, generator, polytopeType, polytopeOptimize, time_ci_alpha, results_ci_alpha, printBody);
        b.run_benchmark();
    } else {
	switch(intersectSet) {
	   case 1: 
	   {
              Benchmark_intersectCoord b("intersectCoord", r, true, warmup, generator, polytopeType, polytopeOptimize, time_ci_alpha, results_ci_alpha, printBody, intersectCoord_intersect, intersectCoord_update);
              b.run_benchmark();
	   } break;
	   case 4:
	   {
              Benchmark_intersectCoord4 b("intersectCoord", r, true, warmup, generator, polytopeType, polytopeOptimize, time_ci_alpha, results_ci_alpha, printBody, intersectCoord_intersect, intersectCoord_update);
              b.run_benchmark();
	   } break;
	   case 8:
	   {
              Benchmark_intersectCoord8 b("intersectCoord", r, true, warmup, generator, polytopeType, polytopeOptimize, time_ci_alpha, results_ci_alpha, printBody, intersectCoord_intersect, intersectCoord_update);
              b.run_benchmark();
	   } break;
	   default:
	      assert(false && "not implemented for this intersectSet!");
	      break;
	}
    }
}
