#include "benchmark.hpp"
#include "../src/volume/volume_helper.hpp"
#include "../polyvest/vol.h"


class Benchmark_App : public Benchmark_base {
    public:
        Benchmark_App(std::string name, int reps, bool convergence, int warmup_reps, const std::string &generator, int polytopeType, const bool polytopeOptimize, const double time_ci_alpha_, const double results_ci_alpha_, const bool printBody)
		: Benchmark_base(name, reps, convergence, warmup_reps, time_ci_alpha_, results_ci_alpha_), generator(generator), polytopeType(polytopeType), polytopeOptimize(polytopeOptimize), printBody(printBody) {}

    protected:
        void initialize () {
            std::cout << "initializing App data..." << std::endl;
        prng_init();
       
       	    solved_body = solved_body_generator()->get(generator,false);
	    if(printBody) {solved_body->print();}

            input = VolumeAppInput_new(solved_body->n, solved_body->bcount);
            
	    input->vol_polytopeType = polytopeType;

	    for(int b=0;b<solved_body->bcount;b++) {
	       input->body[b] = solved_body->body[b];
	       input->type[b] = solved_body->type[b];
	    }
	}
        void reset () {
            // nothing to reset
	}
        double run () {
	    ArbitraryExpNum res = volume_app_ref(input);
		
	    FT exact = solved_body->volume;
	    return (res.num - exact)/exact;
	}
	void finalize() {
	    pc_stack().reset();
            pc_stack().log(0,0,"TODO");
	    //{
            //   PC_Frame<volume_cost_f> frame((void*)volume);
            //   frame.costf()(solved_body->n, solved_body->bcount, (const void**)solved_body->body, (const Body_T**)solved_body->type);
            //}
            pc_stack().print();
	    pc_flops = pc_stack().flops();
	    pc_bytes = pc_stack().bytes();
	}
    private:
	const std::string generator;
	Solved_Body* solved_body;
	FT r0,r1;
	int polytopeType = 0;
	bool polytopeOptimize;
	const bool printBody;
	VolumeAppInput* input; 
};


class Benchmark_Polyvest_App : public Benchmark_base {
    public:
        Benchmark_Polyvest_App(std::string name, int reps, bool convergence, int warmup_reps, const std::string &generator, const bool polytopeOptimize, const double time_ci_alpha_, const double results_ci_alpha_, const bool printBody)
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
	    reset();
	}
        void reset () {
            polyvest_convert(P, Q);
	    Q->check_planes_off = true;
	}
        double run () {
            Q->Preprocess();
            FT res = Q->EstimateVol(step_size); // rout in the step_size also we use
	    FT exact = solved_body->volume;
	    printf("Volume: %.10e\n",res);
	    return (res - exact)/exact;
	}
	void finalize() {
            cout << "stepsize of " << step_size << "\n";
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
        Benchmark_Polyvest_App b("A1_volume", r, true, warmup, generator, polytopeOptimize, time_ci_alpha, results_ci_alpha, printBody);
        b.run_benchmark();
    } else {
        Benchmark_App b("A1_volume", r, true, warmup, generator, polytopeType, polytopeOptimize, time_ci_alpha, results_ci_alpha, printBody);
        b.run_benchmark();
    }
}
