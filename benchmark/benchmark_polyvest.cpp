#include "benchmark.hpp"

#include "../polyvest/vol.h"
#include "../src/volume/volume_helper.hpp"

extern "C" { // must be included C stlye
#include "../src/volume/volume.h"
}

class Polyvest : public Benchmark_base {
    public:
        Polyvest(std::string name, int reps, bool convergence, int warmup_reps) : Benchmark_base(name, reps, convergence, warmup_reps) {}

    protected:
        void initialize () {
            read_polyvest_p("../polyvest/examples/cube_2", &P);
    
            int n = P->n;
            int m = P->m;
    
            Q = new vol::Polyvest_p(m, n);
            polyvest_convert(P, Q);
        }

        void reset () {
            polyvest_convert(P, Q);
        }

        double run () {
            Q->Preprocess();
            Q->EstimateVol(1600);
        }

    private:
        Polytope *P;
        vol::Polyvest_p *Q;
};

int main(int argc, char *argv[]){
    CLI cli(argc,argv,"benchmark");
    CLIFunctionsVolume cliFun(cli);

    int r = 100;
    cliFun.claimOpt('b',"Benchmarking configuration");
    cliFun.add(new CLIF_OptionNumber<int>(&r,'b',"r","100", 1, 100000));
  
    cliFun.preParse();
    if (!cli.parse()) {return -1;}
    cliFun.postParse();

    Polyvest *benchmark = new Polyvest("Polyvest", r, false, 0);
    benchmark->run_benchmark();
}
