#include "benchmark.hpp"

#include "../polyvest/vol.h"
#include "../src/volume/volume_helper.hpp"

extern "C" { // must be included C stlye
#include "../src/volume/preprocess.h"
}

class Polyvest : public Benchmark_base {
    public:
        Polyvest(std::string name, int reps, bool convergence) : Benchmark_base(name, reps, convergence) {}

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

    cli.addOption('r', "100", "number of repetitions");
  
    cliFun.preParse();
    if (!cli.parse()) {return -1;}
    cliFun.postParse();

    int reps = std::stoi(cli.option('r'));

    Polyvest *benchmark = new Polyvest("Polyvest", reps, false);
    benchmark->run_benchmark();
}
