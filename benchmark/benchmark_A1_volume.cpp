#include "benchmark.hpp"
#include "../src/volume/volume_helper.hpp"


class Benchmark_A1 : public Benchmark_base {
    public:
        Benchmark_A1(std::string name, int reps, bool convergence, int n, const std::string &generator) : Benchmark_base(name, reps, convergence), n(n), generator(generator) {}

    protected:
        void initialize () {
            std::cout << "initializing A1 data..." << std::endl;
            
	    if(generator.compare("cube") == 0) {
		bcount = 1;
	        body = (void**)malloc(bcount*sizeof(void*));
	        type = (Body_T**)malloc(bcount*sizeof(Body_T*));
                body[0] = Polytope_new_box(n,1);
		type[0] = &Polytope_T;
		r0 = 1.0;
		r1 = std::sqrt(n);
	    }else if(generator.compare("sphere") == 0) {
		bcount = 1;
	        body = (void**)malloc(bcount*sizeof(void*));
	        type = (Body_T**)malloc(bcount*sizeof(Body_T*));
		FT* center = (FT*)malloc(n*sizeof(FT));
		for(int i=0; i<n; i++) {center[i] = 0;}; center[0] = 1;
                body[0] = Sphere_new(n,2,center);
		type[0] = &Sphere_T;
		r0 = 1.0;
		r1 = 3.0;
	
	    } else {
	        std::cout << "Error: did not find generator " << generator << "\n";
		assert(false);
	    }
        }
        void reset () {
            // nothing to reset
	}
        double run () {
            return volume_ref(n, r0, r1, bcount, (const void**)body, (const Body_T**)type);
	}
    private:
	const std::string generator;
	int n;
	int bcount;
	void** body;
	Body_T** type;
	FT r0,r1;
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
