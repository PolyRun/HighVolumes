#include "benchmark.hpp"
#include "../src/volume/volume_helper.hpp"


class Benchmark_A1 : public Benchmark_base {
    public:
        Benchmark_A1(std::string name, int reps, bool convergence, int warmup_reps, int n, const std::string &generator) : Benchmark_base(name, reps, convergence, warmup_reps), n(n), generator(generator) {}

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
	    }else if(generator.compare("cubeT") == 0) {
		bcount = 1;
	        body = (void**)malloc(bcount*sizeof(void*));
	        type = (Body_T**)malloc(bcount*sizeof(Body_T*));
                body[0] = PolytopeT_new_box(n,1);
		type[0] = &PolytopeT_T;
		r0 = 1.0;
		r1 = std::sqrt(n);
	    }else if(generator.compare("sphere") == 0) {
		bcount = 1;
	        body = (void**)malloc(bcount*sizeof(void*));
	        type = (Body_T**)malloc(bcount*sizeof(Body_T*));
		Ellipsoid* e = Ellipsoid_new(n);
                for(int i=0;i<n;i++) {
                    FT* Ai = Ellipsoid_get_Ai(e,i);
                    Ai[i] = 1.0/4.0;
		    e->a[i] = (i==0);
                }
		body[0] = e;
		type[0] = &Ellipsoid_T;
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
            return volume(n, r0, r1, bcount, (const void**)body, (const Body_T**)type);
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
    
    int n = 20;
    int r = 100;
    cliFun.claimOpt('b',"Benchmarking configuration");
    cliFun.add(new CLIF_OptionNumber<int>(&n,'b',"n","20", 1, 100));
    cliFun.add(new CLIF_OptionNumber<int>(&r,'b',"r","100", 1, 100000));
    
    std::string generator = "cube";
    cliFun.add(new CLIF_Option<std::string>(&generator,'b',"generator","cube", std::map<std::string, std::string>{
                                                     {"cube","cube"},
                                                     {"cubeT","cubeT"},
						     {"sphere","sphere"} }));

    cliFun.preParse();
    if (!cli.parse()) {return -1;}
    cliFun.postParse();

    Benchmark_A1 b("A1_volume", r, true, 0, n, generator);
    b.run_benchmark();
}
