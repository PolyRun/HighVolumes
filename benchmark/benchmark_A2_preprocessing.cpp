#include "benchmark.hpp"
#include "../src/volume/volume_helper.hpp"
#include "../polyvest/vol.h"


class Benchmark_A2 : public Benchmark_base {
    public:
        Benchmark_A2(std::string name, int reps, bool convergence, int warmup_reps, const std::string &generator, int polytopeType, const bool polytopeOptimize, const double time_ci_alpha_, const double results_ci_alpha_, const bool printBody)
		: Benchmark_base(name, reps, convergence, warmup_reps, time_ci_alpha_, results_ci_alpha_), generator(generator), polytopeType(polytopeType), polytopeOptimize(polytopeOptimize), printBody(printBody) {}

    protected:
        void initialize () {
            std::cout << "initializing A2 data..." << std::endl;
        prng_init();
       
       	    solved_body = solved_body_generator()->get(generator,false);
	    if(polytopeOptimize) {
	       solved_body->optimize();
	    }

	    switch(polytopeType) {
            case 1: // column major
                solved_body->polytopeTranspose();
                break;
            case 0: // row major
                break;
            case 2: // CSC format
                solved_body->polytopeCSC();
                break;
            case 3: // JIT format
                solved_body->polytopeJIT();
                break;
	    }
	    if(printBody) {solved_body->print();}
	    //assert(!solved_body->is_normalized);

	    body_out = (void**)malloc(solved_body->bcount * sizeof(void*));
	    for(int c=0;c<solved_body->bcount;c++) {
	       body_out[c] = solved_body->type[c]->clone(solved_body->body[c]);
	    }
	}
        void reset () {
            // nothing to reset
	}
        double run () {
	    ArbitraryExpNum det;
	    preprocess_generic(
			    solved_body->n,
			    solved_body->bcount,
			    (const void**)solved_body->body,
			    body_out,
			    (const Body_T**)solved_body->type,
			    &det
			    );
	    return 0;
	}
	void finalize() {
	    // print:
	    if(printBody) {
               for(int c=0;c<solved_body->bcount;c++) {
	          solved_body->type[c]->print(body_out[c]);
	       }
	    }

	    // pc:
	    pc_stack().reset();
            //{
            //   PC_Frame<volume_cost_f> frame((void*)volume);
            //   frame.costf()(solved_body->n, solved_body->bcount, (const void**)solved_body->body, (const Body_T**)solved_body->type);
            //}
            pc_stack().log(0,0,"TODO");
	    pc_stack().print();
	    pc_flops = pc_stack().flops();
	    pc_bytes = pc_stack().bytes();
	}
    private:
	const std::string generator;
	Solved_Body* solved_body;
	int polytopeType = 0;
	bool polytopeOptimize;
	const bool printBody;
	void ** body_out;
};


class Benchmark_Polyvest_Prep : public Benchmark_base {
    public:
        Benchmark_Polyvest_Prep(std::string name, int reps, bool convergence, int warmup_reps, const std::string &generator, const bool polytopeOptimize, const double time_ci_alpha_, const double results_ci_alpha_, const bool printBody)
		: Benchmark_base(name, reps, convergence, warmup_reps, time_ci_alpha_, results_ci_alpha_), generator(generator), polytopeOptimize(polytopeOptimize), printBody(printBody) {}

    protected:
        void initialize () {
            std::cout << "initializing Polyvest data..." << std::endl;

            solved_body = solved_body_generator()->get(generator,false);
            assert(solved_body->bcount == 1 && "Can maximally have one body for Polyvest.");
            assert(solved_body->type[0] == &Polytope_T && "Can only have polytopes for Polyvest.");

	    if(polytopeOptimize) {
	       solved_body->optimize();
	    }
	    if(printBody) {solved_body->print();}

	    P = (Polytope*)solved_body->body[0];
	     
            int n = P->n;
            int m = P->m;
    
            Q = new vol::Polyvest_p(m, n);
            polyvest_convert(P, Q);
	    Q->check_planes_off = true;
	}
        void reset () {
            delete Q;
            Q = new vol::Polyvest_p(P->m, P->n);
	    polyvest_convert(P, Q);
	    Q->check_planes_off = true;
	}
        double run () {
            Q->Preprocess();
	}
	void finalize() {
	    pc_stack().reset();
            {
               pc_stack().log(0,0,"TODO");
	       //PC_Frame<volume_cost_f> frame((void*)volume);
               //frame.costf()(solved_body->n, solved_body->bcount, (const void**)solved_body->body, (const Body_T**)solved_body->type);
            }
            pc_stack().print();
	    pc_flops = pc_stack().flops();
	    pc_bytes = pc_stack().bytes();
	}
    private:
	const std::string generator;
	Solved_Body* solved_body;
        Polytope* P;
	vol::Polyvest_p *Q;
	bool polytopeOptimize;
	const bool printBody;
};


int main(int argc, char *argv[]){
    CLI cli(argc,argv,"benchmark");
    CLIFunctionsVolume cliFun(cli);
    
    int r = 100;
    int warmup = 0;
    double time_ci_alpha;
    double results_ci_alpha;
    cliFun.claimOpt('b',"Benchmarking configuration");
    cliFun.add(new CLIF_OptionNumber<int>(&r,'b',"r","100", 1, 100000));
    cliFun.add(new CLIF_OptionNumber<int>(&warmup,'b',"warmup","0", 0, 100000));
    cliFun.add(new CLIF_OptionNumber<double>(&time_ci_alpha,'b',"time_ci_alpha","0.95", 0, 1));
    cliFun.add(new CLIF_OptionNumber<double>(&results_ci_alpha,'b',"results_ci_alpha","0.95", 0, 1));

    std::string generator = "cube";
    auto &gen_map = solved_body_generator()->gen_map();
    cliFun.add(new CLIF_Option<std::string>(&generator,'b',"generator","cube_r1.0_10", gen_map));
   
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
                                     {"4",{4, "Polyvest: alternative lib, only for single body polytopes - will preprocess first!"}},
                                    }));

    cliFun.preParse();
    if (!cli.parse()) {return -1;}
    cliFun.postParse();
    
    if(polytopeType==4) {
        Benchmark_Polyvest_Prep b("A1_volume", r, true, warmup, generator, polytopeOptimize, time_ci_alpha, results_ci_alpha, printBody);
        b.run_benchmark();
    } else {
        Benchmark_A2 b("A1_volume", r, true, warmup, generator, polytopeType, polytopeOptimize, time_ci_alpha, results_ci_alpha, printBody);
        b.run_benchmark();
    }
}
