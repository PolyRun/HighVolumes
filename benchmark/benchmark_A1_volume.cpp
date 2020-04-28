#include "benchmark.hpp"
#include "../src/volume/volume_helper.hpp"


class Benchmark_A1 : public Benchmark_base {
    public:
        Benchmark_A1(std::string name, int reps, bool convergence, int warmup_reps, const std::string &generator) : Benchmark_base(name, reps, convergence, warmup_reps), generator(generator) {}

    protected:
        void initialize () {
            std::cout << "initializing A1 data..." << std::endl;
            
            solved_body = solved_body_generator()->get(generator);
	    assert(solved_body->is_preprocessed);
	    r0 = 1.0;
	    r1 = 2*solved_body->n;
	}
        void reset () {
            // nothing to reset
	}
        double run () {
            return volume(solved_body->n, r0, r1, solved_body->bcount, (const void**)solved_body->body, (const Body_T**)solved_body->type)-solved_body->volume;
	}
	void finalize() {
	    pc_stack().reset();
            {
               PC_Frame<volume_cost_f> frame((void*)volume);
               frame.costf()(solved_body->n, solved_body->bcount, (const void**)solved_body->body, (const Body_T**)solved_body->type);
            }
            pc_stack().print();
	    pc_flops = pc_stack().flops();
	    pc_bytes = pc_stack().bytes();
	}
    private:
	const std::string generator;
	Solved_Body* solved_body;
	FT r0,r1;
};

int main(int argc, char *argv[]){
    CLI cli(argc,argv,"benchmark");
    CLIFunctionsVolume cliFun(cli);
    
    int r = 100;
    cliFun.claimOpt('b',"Benchmarking configuration");
    cliFun.add(new CLIF_OptionNumber<int>(&r,'b',"r","100", 1, 100000));
    
    std::string generator = "cube";
    auto &gen_names = solved_body_generator()->identity();
    cliFun.add(new CLIF_Option<std::string>(&generator,'b',"generator","cube_0.5_10", gen_names));

    cliFun.preParse();
    if (!cli.parse()) {return -1;}
    cliFun.postParse();

    Benchmark_A1 b("A1_volume", r, true, 0, generator);
    b.run_benchmark();
}
