#include "benchmark.hpp"
#include "../src/volume/volume_helper.hpp"


class Benchmark_A1 : public Benchmark_base {
    public:
        Benchmark_A1(std::string name, int reps, bool convergence, int n, const std::string &generator) : Benchmark_base(name, reps, convergence), n(n), generator(generator) {}

    protected:
        void initialize () {
            std::cout << "initializing A1 data..." << std::endl;
            
	    if(generator.compare("cube") == 0) {
	        body = (void**)malloc(bcount*sizeof(void*));
	        type = (Body_T**)malloc(bcount*sizeof(Body_T*));

	    } else {
	        std::cout << "Error: did not find generator " << generator << "\n";
		assert(false);
	    }
        }
        void reset () {
            std::cout << "resetting macro benchmark test" << std::endl;
        }
        double run () {
            std::cout << "running macro benchmark test" << std::endl;
        }
    private:
	const std::string generator;
	int n;
	int bcount;
	void** body;
	Body_T** type;
};

int main(int argc, char *argv[]){
    CLI cli(argc,argv,"benchmark");
    CLIFunctionsVolume cliFun(cli);

    cli.addOption('r', "100", "number of repetitions");
    
    int n = 20;
    cliFun.claimOpt('b',"Benchmarking configuration");
    cliFun.add(new CLIF_OptionNumber<int>(&n,'b',"n","20", 1, 100));
    
    std::string generator = "cube";
    cliFun.add(new CLIF_Option<std::string>(&generator,'b',"generator","cube", std::map<std::string, std::string>{
                                                     {"cube","cube"},
						     {"sphere","sphere"} }));

    cliFun.preParse();
    if (!cli.parse()) {return -1;}
    cliFun.postParse();

    int reps = std::stoi(cli.option('r'));

    Benchmark_A1 b("A1_volume", reps, false, n, generator);
    b.run_benchmark();
}
