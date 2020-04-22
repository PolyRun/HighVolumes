#include "benchmark.hpp"
#include "../src/volume/volume_helper.hpp"
#include "../src/util/timer.hpp"

extern "C" {
#include "../src/volume/preprocess.h"
}

#define POLYEXP_BASE ((string) "../../polyvest/examples/")


class Benchmark_A1 : public Benchmark_base {
public:
    Benchmark_A1(int reps, const std::string &polytope_path, Timer_generic *timer) :
        Benchmark_base("A1", reps, false, 0, timer),
        polytope_path(polytope_path)
    {}

protected:
    void initialize () {
        std::cout << "initializing preprocessing:\n";
        std::cout << "read in polytope...\n";

        Polytope *P;
        int err = read_polyvest_p(polytope_path, &P);
        assert(!err &&
               "couldn't read example polytope");

        n = P->n;
        int m = P->m;
        Polytope *Q;
        
        std::cout << "preprocess ellipsoid...\n";
        preprocess(P, &Q, &det);

        bodies = (void **) malloc(sizeof(void *));
        types = (Body_T **) malloc (sizeof(Body_T *));

        bodies[0] = P;
        types[0] = &Polytope_T;
    }
    
    void reset () {}

    double run () {
        return volume_ref(n, 1, 2*n, 1, (const void **) bodies, (const Body_T **) types) * det;
    }
    
private:
    std::string polytope_path;
    void** bodies;
    Body_T** types;
    int n;
    FT det;
};


int main(int argc, char *argv[]){
    CLI cli(argc,argv,"benchmark");
    CLIFunctionsVolume cliFun(cli);

    cli.addOption('r', "5", "number of repetitions");
    
    cli.addOption('P', "cube_2", "input polytope");

    cli.addOption('t', "1", "timer: 0 chrono, 1 tsc");
    
    cliFun.preParse();
    if (!cli.parse()) {return -1;}
    cliFun.postParse();
    

    std::string path_from_exec = cli.getPathFromExec() + POLYEXP_BASE + cli.option('P');

    int reps = std::stoi(cli.option('r'));

    Timer_generic *timer;
    if (std::stoi(cli.option('t')) == 0) {
        Timer t = Timer();
        timer = (Timer_generic *) &t;
    }
    else {
        Tsc t = Tsc();
        timer = (Timer_generic *) &t;
    }
    //Tsc clocks = Tsc();
    Benchmark_A1 b(reps, path_from_exec, timer);
    b.run_benchmark();
}
