#include "benchmark.hpp"
#include "../src/volume/volume_helper.hpp"


class Benchmark_test : public Benchmark_base {
    public:
        Benchmark_test(std::string name, int reps, bool convergence, int warmup_reps, const double time_ci_alpha_, const double results_ci_alpha_, int n)
		: Benchmark_base(name, reps, convergence, warmup_reps, time_ci_alpha_, results_ci_alpha_), n(n) {}

    protected:
        void initialize () {
            std::cout << "initializing test data..." << std::endl;
	    
            jit_clear();
            func = (double (*)(double*,double*)) jit_head();
	    jit_Table_32* t32 = NULL;// empty list
	    
	    t32 = jit_immediate_32_via_data(0,0,0,0, 0, t32);
	    t32 = jit_immediate_32_via_data(0,0,0,0, 4, t32);
	    
	    for(int i=0;i<n;i++) {
	       //int ii = prng_get_random_int_in_range(-10,10);
               //ii = std::min(4*n-4, std::max(0, i*4 + ii));
	       //t32 = jit_immediate_32_via_data(2.0,2.0,2.0,2.0, 2, t32);
	       jit_loadu_ymm(jit_rsi,4*8*i,2);
	       jit_vmulpd_mem_ymm(jit_rdi,4*8*i,2,1);
               jit_vmaxpd_ymm(1,4*(i%2),4*(i%2));
               //jit_vmaxpd_ymm(1,0,0);
	    }
            
	    jit_vmaxpd_ymm(4,0,0);
	    jit_emit_vzeroupper();
	    
	    jit_emit_return();
            bytes_op = jit_head() - (uint8_t*)func;
	    
	    jit_table_32_consume(t32);
            
	    bytes_data = jit_head() - (uint8_t*)func - bytes_op;
	    
	    //jit_print();
     
	    x = (FT*)(aligned_alloc(32, 4*n*sizeof(FT))); // align this to 32
	    y = (FT*)(aligned_alloc(32, 4*n*sizeof(FT))); // align this to 32
	    for(int i=0;i<4*n;i++) {x[i] = i; y[i]=2;}

            double res = func(x,y);
	    std::cout << "res: " << res << " " << 4*n << "\n";

	    reset();
        }
        void reset () {
	}
        double run () {
            double res = func(x,y);
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
	double* y;
	double (*func)(double*,double*);
	size_t bytes_op;
	size_t bytes_data;
};

class Benchmark_testa : public Benchmark_base {
    public:
        Benchmark_testa(std::string name, int reps, bool convergence, int warmup_reps, const double time_ci_alpha_, const double results_ci_alpha_, int n)
		: Benchmark_base(name, reps, convergence, warmup_reps, time_ci_alpha_, results_ci_alpha_), n(n) {}

    protected:
        void initialize () {
            std::cout << "initializing test data..." << std::endl;
	    
            jit_clear();
            func = (double (*)(double*,double*)) jit_head();
	    jit_Table_8* t8 = NULL;// empty list
	    
	    t8 = jit_immediate_8_via_data(0, 0, t8);
	    t8 = jit_immediate_8_via_data(0, 4, t8);
	    
	    for(int i=0;i<4*n;i++) {
	       jit_load_sd(jit_rsi,8*i,2);
	       jit_vmulsd_mem(jit_rdi,8*i,2,1);
               jit_vmaxsd(1,4*(i%2),4*(i%2));
	    }
            
	    jit_vmaxsd(4,0,0);
	    //jit_emit_vzeroupper();
	    
	    jit_emit_return();
            bytes_op = jit_head() - (uint8_t*)func;
	    
	    jit_table_8_consume(t8);
            
	    bytes_data = jit_head() - (uint8_t*)func - bytes_op;
	    
	    //jit_print();
     
	    x = (FT*)(aligned_alloc(32, 4*n*sizeof(FT))); // align this to 32
	    y = (FT*)(aligned_alloc(32, 4*n*sizeof(FT))); // align this to 32
	    for(int i=0;i<4*n;i++) {x[i] = i; y[i]=2;}

            double res = func(x,y);
	    std::cout << "res: " << res << " " << 4*n << "\n";

	    reset();
        }
        void reset () {
	}
        double run () {
            double res = func(x,y);
	    //std::cout << "res: " << res << "\n";
	    return 0;
	}
	void finalize() {
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
	double* y;
	double (*func)(double*,double*);
	size_t bytes_op;
	size_t bytes_data;
};



class Benchmark_test2 : public Benchmark_base {
    public:
        Benchmark_test2(std::string name, int reps, bool convergence, int warmup_reps, const double time_ci_alpha_, const double results_ci_alpha_, int n)
		: Benchmark_base(name, reps, convergence, warmup_reps, time_ci_alpha_, results_ci_alpha_), n(n) {}

    protected:
        void initialize () {
            std::cout << "initializing test data..." << std::endl;
	    
            jit_clear();
            func = (double (*)(double*,double*)) jit_head();
	    jit_Table_32* t32 = NULL;// empty list
	    
	    t32 = jit_immediate_32_via_data(0,0,0,0, 0, t32);
	    t32 = jit_immediate_32_via_data(0,0,0,0, 4, t32);
	    
	    for(int i=0;i<n;i++) {
	       //int ii = prng_get_random_int_in_range(-10,10);
               //ii = std::min(4*n-4, std::max(0, i*4 + ii));
	       //t32 = jit_immediate_32_via_data(2.0,2.0,2.0,2.0, 2, t32);
	       jit_loadu_ymm(jit_rsi,4*8*i,2);
	       int ii = prng_get_random_int_in_range(0,4*n-4);
	       jit_vmulpd_mem_ymm(jit_rdi,ii*8,2,1);
               jit_vmaxpd_ymm(1,4*(i%2),4*(i%2));
               //jit_vmaxpd_ymm(1,0,0);
	    }
            
	    jit_vmaxpd_ymm(4,0,0);
	    jit_emit_vzeroupper();
	    
	    jit_emit_return();
            bytes_op = jit_head() - (uint8_t*)func;
	    
	    jit_table_32_consume(t32);
            
	    bytes_data = jit_head() - (uint8_t*)func - bytes_op;
	    
	    //jit_print();
     
	    x = (FT*)(aligned_alloc(32, 4*n*sizeof(FT))); // align this to 32
	    y = (FT*)(aligned_alloc(32, 4*n*sizeof(FT))); // align this to 32
	    for(int i=0;i<4*n;i++) {x[i] = i; y[i]=2;}

            double res = func(x,y);
	    std::cout << "res: " << res << " " << 4*n << "\n";

	    reset();
        }
        void reset () {
	}
        double run () {
            double res = func(x,y);
	    //std::cout << "res: " << res << "\n";
	    return 0;
	}
	void finalize() {
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
	double* y;
	double (*func)(double*,double*);
	size_t bytes_op;
	size_t bytes_data;
};



class Benchmark_test2a : public Benchmark_base {
    public:
        Benchmark_test2a(std::string name, int reps, bool convergence, int warmup_reps, const double time_ci_alpha_, const double results_ci_alpha_, int n)
		: Benchmark_base(name, reps, convergence, warmup_reps, time_ci_alpha_, results_ci_alpha_), n(n) {}

    protected:
        void initialize () {
            std::cout << "initializing test data..." << std::endl;
	    
            jit_clear();
            func = (double (*)(double*,double*)) jit_head();
	    jit_Table_8* t8 = NULL;// empty list
	    
	    t8 = jit_immediate_8_via_data(0, 0, t8);
	    t8 = jit_immediate_8_via_data(0, 4, t8);
	    
	    for(int i=0;i<4*n;i++) {
	       jit_load_sd(jit_rsi,8*i,2);
	       int ii = prng_get_random_int_in_range(0,4*n-1);
	       jit_vmulsd_mem(jit_rdi,ii*8,2,1);
               jit_vmaxsd(1,4*(i%2),4*(i%2));
	    }
            
	    jit_vmaxsd(4,0,0);
	    //jit_emit_vzeroupper();
	    
	    jit_emit_return();
            bytes_op = jit_head() - (uint8_t*)func;
	    
	    jit_table_8_consume(t8);
            
	    bytes_data = jit_head() - (uint8_t*)func - bytes_op;
	    
	    //jit_print();
     
	    x = (FT*)(aligned_alloc(32, 4*n*sizeof(FT))); // align this to 32
	    y = (FT*)(aligned_alloc(32, 4*n*sizeof(FT))); // align this to 32
	    for(int i=0;i<4*n;i++) {x[i] = i; y[i]=2;}

            double res = func(x,y);
	    std::cout << "res: " << res << " " << 4*n << "\n";

	    reset();
        }
        void reset () {
	}
        double run () {
            double res = func(x,y);
	    //std::cout << "res: " << res << "\n";
	    return 0;
	}
	void finalize() {
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
	double* y;
	double (*func)(double*,double*);
	size_t bytes_op;
	size_t bytes_data;
};

class Benchmark_test3 : public Benchmark_base {
    public:
	typedef double (*FF)(const double*);
        Benchmark_test3(std::string name, int reps, bool convergence, int warmup_reps, const double time_ci_alpha_, const double results_ci_alpha_, int n, int m, int w)
		: Benchmark_base(name, reps, convergence, warmup_reps, time_ci_alpha_, results_ci_alpha_), n(n), m(m),w(w) {}

    protected:
        void initialize () {
            std::cout << "initializing test data..." << std::endl;
	    
            jit_clear();
            funcs.resize(m,NULL);
	    jit_Table_32* t32 = NULL;// empty list
	    for(int j=0;j<m;j++) {
	       funcs[j] = (double (*)(const double*)) jit_head();
	       t32 = jit_immediate_32_via_data(0,0,0,0, 0, t32);
	       t32 = jit_immediate_32_via_data(0,0,0,0, 4, t32);
	    
 	       for(int i=0;i<w;i++) {
	          t32 = jit_immediate_32_via_data(j,j,j,j, 2, t32);
	          int ii = prng_get_random_int_in_range(0,4*n-4);
	          jit_vmulpd_mem_ymm(jit_rdi,ii*8,2,1);
                  jit_vmaxpd_ymm(1,4*(i%2),4*(i%2));
	       }
	       jit_vmaxpd_ymm(4,0,0);
	       jit_emit_vzeroupper();
	       jit_emit_return();
	    }
            bytes_op = jit_head() - (uint8_t*)funcs[0];
	    jit_table_32_consume(t32);
	    bytes_data = jit_head() - (uint8_t*)funcs[0] - bytes_op;
	    
	    //jit_print();
     
	    x = (FT*)(aligned_alloc(32, 4*n*sizeof(FT))); // align this to 32
	    for(int i=0;i<4*n;i++) {x[i] = i;}

            double res = funcs[0](x);
	    std::cout << "res: " << res << " " << 4*n << "\n";

	    reset();
        }
        void reset () {
	   d = prng_get_random_int_in_range(0,m-1);
	}
        double run () {
            double res = funcs[d](x);
	    //std::cout << "res: " << res << "\n";
	    return 0;
	}
	void finalize() {
	    pc_stack().reset();
            {
                pc_stack().log(8*w, 8*w*sizeof(FT)," 4 mul, 4 max, ld x and consts!");
	    }
            pc_stack().print();
	    pc_flops = pc_stack().flops();
	    pc_bytes = pc_stack().bytes();
	}
    protected:
        int n;
        int m;
        int w;
	double* x;
	std::vector<FF> funcs;
	//double (*func)(double*,double*);
	size_t bytes_op;
	size_t bytes_data;
	int d=0;// dimension
};


class Benchmark_test4 : public Benchmark_base {
    public:
	typedef double (*FF)(const double*);
        Benchmark_test4(std::string name, int reps, bool convergence, int warmup_reps, const double time_ci_alpha_, const double results_ci_alpha_, int n, int m, int w)
		: Benchmark_base(name, reps, convergence, warmup_reps, time_ci_alpha_, results_ci_alpha_), n(n), m(m),w(w) {}

    protected:
        void initialize () {
            std::cout << "initializing test data..." << std::endl;
	    
            jit_clear();
            funcs.resize(m,NULL);
	    jit_Table_16* t16 = NULL;// empty list
	    for(int j=0;j<m;j++) {
	       funcs[j] = (double (*)(const double*)) jit_head();
	       t16 = jit_immediate_16_via_data(0,0, 0, t16);
	       t16 = jit_immediate_16_via_data(0,0, 4, t16);
	    
 	       for(int i=0;i<w;i++) {
	          t16 = jit_immediate_16_via_data(j,j, 2, t16);
	          int ii = prng_get_random_int_in_range(0,2*n-2);
	          jit_vmulpd_mem_xmm(jit_rdi,ii*8,2,1);
                  jit_vmaxpd_xmm(1,4*(i%2),4*(i%2));
	       }
	       jit_vmaxpd_xmm(4,0,0);
	       jit_emit_return();
	    }
            bytes_op = jit_head() - (uint8_t*)funcs[0];
	    jit_table_16_consume(t16);
	    bytes_data = jit_head() - (uint8_t*)funcs[0] - bytes_op;
	    
	    //jit_print();
     
	    x = (FT*)(aligned_alloc(32, 2*n*sizeof(FT))); // align this to 32
	    for(int i=0;i<2*n;i++) {x[i] = i;}

            double res = funcs[0](x);
	    std::cout << "res: " << res << " " << 4*n << "\n";

	    reset();
        }
        void reset () {
	   d = prng_get_random_int_in_range(0,m-1);
	}
        double run () {
            double res = funcs[d](x);
	    //std::cout << "res: " << res << "\n";
	    return 0;
	}
	void finalize() {
	    pc_stack().reset();
            {
                pc_stack().log(4*w, 4*w*sizeof(FT)," 4 mul, 4 max, ld x and consts!");
	    }
            pc_stack().print();
	    pc_flops = pc_stack().flops();
	    pc_bytes = pc_stack().bytes();
	}
    protected:
        int n;
        int m;
        int w;
	double* x;
	std::vector<FF> funcs;
	//double (*func)(double*,double*);
	size_t bytes_op;
	size_t bytes_data;
	int d=0;// dimension
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
    int m = 1;
    int w = 1;
    cliFun.add(new CLIF_OptionNumber<int>(&n,'b',"n","1000", 1, 10000000));
    cliFun.add(new CLIF_OptionNumber<int>(&m,'b',"m","10", 1, 10000000));
    cliFun.add(new CLIF_OptionNumber<int>(&w,'b',"w","100", 1, 10000000));

    std::string experiment = "";
    cliFun.add(new CLIF_Option<std::string>(&experiment,'b',"experiment","test",
                                    {
                                     {"test",{"test", "linear access - n=data ymm"}},
                                     {"testa",{"test", "linear access - n=data sd"}},
                                     {"test2",{"test2", "random access - n=data ymm"}},
                                     {"test2a",{"test2", "random access - n=data sd"}},
                                     {"test3",{"test3", "ymm, many functions, random access - n=data,m=Nfunctions,w=dataPerFunction"}},
                                     {"test4",{"test4", "xmm, many functions, random access - n=data,m=Nfunctions,w=dataPerFunction"}},
                                    }));

    cliFun.preParse();
    if (!cli.parse()) {return -1;}
    cliFun.postParse();
    
    if(experiment == "test") {
        Benchmark_test b(experiment, r, true, warmup, time_ci_alpha, results_ci_alpha, n);
        b.run_benchmark();
    } else if(experiment == "testa") {
        Benchmark_testa b(experiment, r, true, warmup, time_ci_alpha, results_ci_alpha, n);
        b.run_benchmark();
    } else if(experiment == "test2") {
        Benchmark_test2 b(experiment, r, true, warmup, time_ci_alpha, results_ci_alpha, n);
        b.run_benchmark();
    } else if(experiment == "test2a") {
        Benchmark_test2a b(experiment, r, true, warmup, time_ci_alpha, results_ci_alpha, n);
        b.run_benchmark();
    } else if(experiment == "test3") {
        Benchmark_test3 b(experiment, r, true, warmup, time_ci_alpha, results_ci_alpha, n,m,w);
        b.run_benchmark();
    } else if(experiment == "test4") {
        Benchmark_test4 b(experiment, r, true, warmup, time_ci_alpha, results_ci_alpha, n,m,w);
        b.run_benchmark();
    } else {
        assert(false && "experiment not hantled!");
    }
}
