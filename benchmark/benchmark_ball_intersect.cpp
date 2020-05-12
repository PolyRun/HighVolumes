#include "benchmark.hpp"
#include "../src/volume/volume_helper.hpp"


class Benchmark_intersect : public Benchmark_base {
    public:
        Benchmark_intersect(std::string name, int reps, bool convergence, int warmup_reps, int n, const double time_ci_alpha_, const double results_ci_alpha_)
		: Benchmark_base(name, reps, convergence, warmup_reps, time_ci_alpha_, results_ci_alpha_), n(n) {}

    protected:
        void initialize () {
            std::cout << "initializing intersect data..." << std::endl;
	    
	    x = (FT*)(aligned_alloc(32, n*sizeof(FT))); // align this to 32
	    d = (FT*)(aligned_alloc(32, n*sizeof(FT))); // align this to 32
            
	    reset();
        }
        void reset () {
            r = prng_get_random_double_in_range(0.1,10);
	    for(int i=0; i<n;i++) {
	        x[i] = prng_get_random_double_0_1()*r*0.9;
	        d[i] = prng_get_random_double_normal();
	    }
	    dd = prng_get_random_int_in_range(0,n-1);
	}
        double run () {
	    FT t0, t1;
	    Ball_intersect(n, r, x, d, &t0, &t1);
            return t0-t1;
	}
	void finalize() {
	    pc_stack().reset();
            {
                PC_Frame<Ball_intersect_cost_f> frame((void*)Ball_intersect);
                frame.costf()(n);
            }
            pc_stack().print();
	    pc_flops = pc_stack().flops();
	    pc_bytes = pc_stack().bytes();
	}
    protected:
	int n;
	FT* x;
	FT* d;
	FT r;
	int dd;
};

class Benchmark_intersectCoord : public Benchmark_intersect {
    public:
        Benchmark_intersectCoord(std::string name, int reps, bool convergence, int warmup_reps, int n, const double time_ci_alpha_, const double results_ci_alpha_)
		: Benchmark_intersect(name, reps, convergence, warmup_reps, n, time_ci_alpha_, results_ci_alpha_) {}
    
	void finalize() {
	    pc_stack().reset();
            {
                PC_Frame<Ball_intersectCoord_cost_f> frame((void*)Ball_intersectCoord);
                frame.costf()(n);
            }
            pc_stack().print();
	    pc_flops = pc_stack().flops();
	    pc_bytes = pc_stack().bytes();
	}
    protected:
        double run () {
	    FT t0, t1;
	    Ball_intersectCoord(n, r, x, dd, &t0, &t1);
            return 0;
	}
};

int main(int argc, char *argv[]){
    CLI cli(argc,argv,"benchmark");
    CLIFunctionsVolume cliFun(cli);
    
    int n = 20;
    int r = 100;
    int warmup = 0;
    double time_ci_alpha;
    double results_ci_alpha;
    cliFun.claimOpt('b',"Benchmarking configuration");
    cliFun.add(new CLIF_OptionNumber<int>(&n,'b',"n","20", 1, 200));
    cliFun.add(new CLIF_OptionNumber<int>(&r,'b',"r","100", 1, 10000000));
    cliFun.add(new CLIF_OptionNumber<int>(&warmup,'b',"warmup","0", 0, 10000000));
    cliFun.add(new CLIF_OptionNumber<double>(&time_ci_alpha,'b',"time_ci_alpha","0.95", 0, 1));
    cliFun.add(new CLIF_OptionNumber<double>(&results_ci_alpha,'b',"results_ci_alpha","0.95", 0, 1));
 
    std::string intersect = "intersect";
    cliFun.add(new CLIF_Option<std::string>(&intersect,'b',"intersect","intersect", {
                                                     {"intersect",      {"intersect",     "random direction intersection"}},
						     {"intersectCoord", {"intersectCoord","coordinate direction intersection"}} }));

    cliFun.preParse();
    if (!cli.parse()) {return -1;}
    cliFun.postParse();
    
    if(intersect.compare("intersect")==0) {
        Benchmark_intersect b("intersect", r, true, warmup, n, time_ci_alpha, results_ci_alpha);
        b.run_benchmark();
    } else {
        Benchmark_intersectCoord b("intersectCoord", r, true, warmup, n, time_ci_alpha, results_ci_alpha);
        b.run_benchmark();
    }
}
