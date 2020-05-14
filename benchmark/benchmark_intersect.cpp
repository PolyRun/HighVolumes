#include "benchmark.hpp"
#include "../src/volume/volume_helper.hpp"


class Benchmark_intersect : public Benchmark_base {
    public:
        Benchmark_intersect(std::string name, int reps, bool convergence, int warmup_reps, const std::string &generator, const int polytopeType, const bool polytopeOptimize, const double time_ci_alpha_, const double results_ci_alpha_)
		: Benchmark_base(name, reps, convergence, warmup_reps, time_ci_alpha_, results_ci_alpha_), generator(generator), polytopeType(polytopeType), polytopeOptimize(polytopeOptimize) {}

    protected:
        void initialize () {
            std::cout << "initializing intersect data..." << std::endl;
	   
            solved_body = solved_body_generator()->get(generator,false);
	    if(polytopeOptimize) {
	       solved_body->optimize();
	    }

	    switch(polytopeType) {
            case 0: // column major
                solved_body->polytopeTranspose();
                break;
            case 1: // row major
                break;
            case 2: // CSC format
                solved_body->polytopeCSC();
                break;
            case 3: // JIT format
                solved_body->polytopeJIT();
                break;
	    }
	    assert(solved_body->is_normalized);
	    assert(solved_body->bcount == 1);
	    
	    int n = solved_body->n;
	    x = (FT*)(aligned_alloc(32, n*sizeof(FT))); // align this to 32
	    d = (FT*)(aligned_alloc(32, n*sizeof(FT))); // align this to 32
	    
	    int cache_size = solved_body->type[0]->cacheAlloc(solved_body->body[0]);
            cache = aligned_alloc(32, cache_size); // align this to 32

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
                pc_stack().log(0,0, "random double - TODO");
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
	bool polytopeOptimize = false;
};

class Benchmark_intersectCoord : public Benchmark_intersect {
    public:
        Benchmark_intersectCoord(std::string name, int reps, bool convergence, int warmup_reps, const std::string &generator, const int polytopeType, const bool polytopeOptimize, const double time_ci_alpha_, const double results_ci_alpha_)
		: Benchmark_intersect(name, reps, convergence, warmup_reps, generator, polytopeType, polytopeOptimize, time_ci_alpha_, results_ci_alpha_) {}
    
    	void finalize() {
	    pc_stack().reset();
            {
		{
		    PC_Frame<intersect_cost_f> frame((void*)solved_body->type[0]->intersectCoord);
                    frame.costf()(solved_body->body[0]);
		}

                pc_stack().log(0, 0, "random double - TODO");
	        // Reading and writing x[dd] with one add in between
                pc_stack().log(1, 2, "x[dd] += t;");
                
                // body intersectCoord
		{
	            PC_Frame<cacheUpdateCoord_cost_f> frame((void*) solved_body->type[0]->cacheUpdateCoord);
                    frame.costf()(solved_body->body[0]);
		}
	    }
            pc_stack().print();
	    pc_flops = pc_stack().flops();
	    pc_bytes = pc_stack().bytes();
	}
    protected:
        double run () {
	    FT t0, t1;
	    solved_body->type[0]->intersectCoord(solved_body->body[0], x, dd, &t0, &t1, cache);
            
	    // step now
	    FT t = prng_get_random_double_in_range(t0,t1);
            x[dd] += t;
            solved_body->type[0]->cacheUpdateCoord(solved_body->body[0], dd, t, cache);
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
    
    std::string intersect = "intersect";
    cliFun.add(new CLIF_Option<std::string>(&intersect,'b',"intersect","intersect", {
                                                     {"intersect",      {"intersect",     "random direction intersection"}},
						     {"intersectCoord", {"intersectCoord","coordinate direction intersection"}} }));
   
    bool polytopeOptimize = false;
    cliFun.add(new CLIF_Option<bool>(&polytopeOptimize,'b',"polytopeOptimize","false", {
                                                     {"false",{false,"-"}},
						     {"true", {true, "Sort constraints to optimize access pattern"}} }));

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
    
    if(intersect.compare("intersect")==0) {
        Benchmark_intersect b("intersect", r, true, warmup, generator, polytopeType, polytopeOptimize, time_ci_alpha, results_ci_alpha);
        b.run_benchmark();
    } else {
        Benchmark_intersectCoord b("intersectCoord", r, true, warmup, generator, polytopeType, polytopeOptimize, time_ci_alpha, results_ci_alpha);
        b.run_benchmark();
    }
}
