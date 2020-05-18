#include "benchmark.hpp"
#include "../src/volume/volume_helper.hpp"


class Benchmark_test : public Benchmark_base {
    public:
        Benchmark_test(std::string name, int reps, bool convergence, int warmup_reps, const double time_ci_alpha_, const double results_ci_alpha_, int n)
		: Benchmark_base(name, reps, convergence, warmup_reps, time_ci_alpha_, results_ci_alpha_), n(4*n) {}

    protected:
        void initialize () {
            std::cout << "initializing test data..." << std::endl;
	    
            jit_clear();
            func = (double (*)(double*)) jit_head();
	    jit_Table_32* t32 = NULL;// empty list
	    
	    t32 = jit_immediate_32_via_data(0,0,0,0, 0, t32);
	    
	    for(int i=0;i<n;i++) {
	       t32 = jit_immediate_32_via_data(2.0,2.0,2.0,2.0, 1, t32);
	       jit_vmulpd_mem_ymm(jit_rdi,32*i,1,1);
               jit_vmaxpd_ymm(1,0,0);
	    }
	    
	    jit_emit_return();
            bytes_op = jit_head() - (uint8_t*)func;
	    
	    jit_table_32_consume(t32);
            
	    bytes_data = jit_head() - (uint8_t*)func - bytes_op;
	    
	    //jit_print();
     
	    x = (FT*)(aligned_alloc(32, n*sizeof(FT))); // align this to 32
	    for(int i=0;i<n;i++) {x[i] = i;}

	    reset();
        }
        void reset () {
	}
        double run () {
            double res = func(x);
	    //std::cout << "res: " << res << "\n";
	    return 0;
	}
	void finalize() {
	    std::cout << "Bytes  op: " << bytes_op << ", data: "<< bytes_data <<"\n";
	    std::cout << "Expected flops/cycle: mul and max per 3 cycles, times 4 = 8/3 = 2.66\n";
	    std::cout << "Expected bytes/cycle: data in x, but also constants!\n";
	    pc_stack().reset();
            {
                pc_stack().log(8*n, 8*n*sizeof(FT)," 4 mul, 4 max, ld x and consts!");
	    }
            pc_stack().print();
	    pc_flops = pc_stack().flops();
	    pc_bytes = pc_stack().bytes();
	}
    protected:
        int n;
	double* x;
	double (*func)(double*);
	size_t bytes_op;
	size_t bytes_data;
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
    
    int n = 1;
    cliFun.add(new CLIF_OptionNumber<int>(&n,'b',"n","1000", 1, 10000000));

    std::string experiment = "";
    cliFun.add(new CLIF_Option<std::string>(&experiment,'b',"experiment","test",
                                    {
                                     {"test",{"test", "first experiment"}},
                                    }));

    cliFun.preParse();
    if (!cli.parse()) {return -1;}
    cliFun.postParse();
    
    if(experiment == "test") {
        Benchmark_test b(experiment, r, true, warmup, time_ci_alpha, results_ci_alpha, n);
        b.run_benchmark();
    } else {
        assert(false && "experiment not hantled!");
    }
}
