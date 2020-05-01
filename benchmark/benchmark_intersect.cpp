#include "benchmark.hpp"
#include "../src/volume/volume_helper.hpp"


class Benchmark_intersect : public Benchmark_base {
    public:
        Benchmark_intersect(std::string name, int reps, bool convergence, int warmup_reps, const std::string &generator, const bool polytopeTranspose)
		: Benchmark_base(name, reps, convergence, warmup_reps), generator(generator), polytopeTranspose(polytopeTranspose){}

    protected:
        void initialize () {
            std::cout << "initializing intersect data..." << std::endl;
	    
	    solved_body = solved_body_generator()->get(generator, polytopeTranspose);
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
	    FT t0, t1;
	    solved_body->type[0]->intersect(solved_body->body[0], x, d, &t0, &t1);
            return t0-t1;
	}
	void finalize() {
	    pc_stack().reset();
            {
                PC_Frame<intersect_cost_f> frame((void*)solved_body->type[0]->intersect);
                frame.costf()(solved_body->body[0]);
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
	bool polytopeTranspose = false;
};

class Benchmark_intersectCoord : public Benchmark_intersect {
    public:
        Benchmark_intersectCoord(std::string name, int reps, bool convergence, int warmup_reps, const std::string &generator, const bool polytopeTranspose)
		: Benchmark_intersect(name, reps, convergence, warmup_reps, generator, polytopeTranspose) {}
    
    	void finalize() {
	    pc_stack().reset();
            {
                PC_Frame<intersect_cost_f> frame((void*)solved_body->type[0]->intersectCoord);
                frame.costf()(solved_body->body[0]);
            }
            pc_stack().print();
	    pc_flops = pc_stack().flops();
	    pc_bytes = pc_stack().bytes();
	}
    protected:
        double run () {
	    FT t0, t1;
	    solved_body->type[0]->intersectCoord(solved_body->body[0], x, dd, &t0, &t1, cache);
            return 0;
	}
};

int main(int argc, char *argv[]){
    CLI cli(argc,argv,"benchmark");
    CLIFunctionsVolume cliFun(cli);
    
    int r = 100;
    int warmup = 0;
    cliFun.claimOpt('b',"Benchmarking configuration");
    cliFun.add(new CLIF_OptionNumber<int>(&r,'b',"r","100", 1, 10000000));
    cliFun.add(new CLIF_OptionNumber<int>(&warmup,'b',"warmup","0", 0, 10000000));
    
    std::string generator = "cube";
    auto &gen_map = solved_body_generator()->gen_map();
    cliFun.add(new CLIF_Option<std::string>(&generator,'b',"generator","cube_r1.0_10", gen_map));
    
    std::string intersect = "intersect";
    cliFun.add(new CLIF_Option<std::string>(&intersect,'b',"intersect","intersect", {
                                                     {"intersect",      {"intersect",     "random direction intersection"}},
						     {"intersectCoord", {"intersectCoord","coordinate direction intersection"}} }));
    
    bool polytopeTranspose = false;
    cliFun.add(new CLIF_Option<bool>(&polytopeTranspose,'b',"polytopeTranspose","false", {
                                                     {"false",{false, "Polytope format / rows"}},
						     {"true",{true, "PolytopeT format / columns"}} }));

    cliFun.preParse();
    if (!cli.parse()) {return -1;}
    cliFun.postParse();
    
    if(intersect.compare("intersect")==0) {
        Benchmark_intersect b("intersect", r, true, warmup, generator, polytopeTranspose);
        b.run_benchmark();
    } else {
        Benchmark_intersectCoord b("intersectCoord", r, true, warmup, generator, polytopeTranspose);
        b.run_benchmark();
    }
}
