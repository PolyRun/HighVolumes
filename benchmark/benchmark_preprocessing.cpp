#include "benchmark.hpp"
#include "../src/volume/volume_helper.hpp"
#include "../src/util/timer.hpp"

extern "C" {
#include "../src/volume/volume.h"
}

#define POLYEXP_BASE ((string) "../../polyvest/examples/")


class Benchmark_preprocessing : public Benchmark_base {
public:
    Benchmark_preprocessing(std::string name,
                            int reps,
                            bool convergence,
                            int warmup_reps,
                            const std::string &polytope_path,
                            Timer_generic *timer) :
        Benchmark_base(name, reps, convergence, warmup_reps, 0.95,0.95, timer),
        polytope_path(polytope_path)
    {}

protected:
    void initialize () {
        std::cout << "initializing preprocessing:\n";
        std::cout << "read in polytope...\n";

        int err = read_polyvest_p(polytope_path, &P);
        assert(!err &&
               "couldn't read example polytope");

        int n = P->n;
        int m = P->m;
        
        std::cout << "initialize ellipsoid...\n";
        init_ellipsoid(P, &R2, &ori);

        std::cout << "R2: " << R2 << "\nori: ";
        for (int i = 0; i < n; i++) std::cout << ori[i] << " ";
        std::cout << std::endl;
            
        std::cout << "get opcount...\n";
        preprocess_opcount(P, R2, ori, &Q, &det, &c1, &c2, &c3, &c4);

        
        nadds = c1 * (2*n*n + n) +
            c2 * (n+1) +
            c3 * (n*n+n+1) +
            c4 * (n*n+1) +
            (n*n*n + 6*n*n + 5*n)/6 + m*n*(n+1)/2 + 7;

        nmults = c1 * (4*n*n + n) +
            c2 * n +
            c3 * (n*n+n+2) +
            c4 * (n*n + n) +
            (n*n*n + 3*n*n - 4*n)/6 + m*n*(n+1)/2 + 2*n + m + 9;

        ndivs = c1 * n + (n*n + n)/2 + 8;

        nsqrts = c1 * n + n;
        

        std::cout << "number of cuts: " << c1
                  << "\niterations of first inner loop: " << c2
                  << "\niterations of second inner loop: " << c3
                  << "\nnadds: " << nadds
                  << "\nnmults: " << nmults
                  << "\nndivs: " << ndivs
                  << "\nsqrts: " << nsqrts << std::endl;
    }
    
    void reset () {}

    double run () {
        preprocess(P, &Q, &det);
        return 0;
    }
    
private:
    // input
    Polytope *P; // the input polytope
    FT R2; // radius squared of initial ellipsoid
    FT *ori; // center of initiali ellipsoid

    // output
    Polytope *Q; // the output polytope
    FT det; // the determinant of the transformation Q -> P
    int c1, c2, c3, c4; // variables holding counts of loop iterations in preprocess on this input
    int nmults, nadds, ndivs, nsqrts; // total opcount of preprocess on this input
    std::string polytope_path;
};

int main(int argc, char *argv[]){
    CLI cli(argc,argv,"benchmark");
    CLIFunctionsVolume cliFun(cli);

    int r = 100;
    cliFun.claimOpt('b',"Benchmarking configuration");
    cliFun.add(new CLIF_OptionNumber<int>(&r,'b',"r","100", 1, 100000));
    
    cli.addOption('P', "cube_2", "input polytope");

    cli.addOption('t', "1", "timer: 0 chrono, 1 tsc");
    
    cliFun.preParse();
    if (!cli.parse()) {return -1;}
    cliFun.postParse();
    
    std::string path = cli.getPath();
    std::string path_from_exec = "";
    reverse(path.begin(), path.end());
    size_t pos = path.find('/');
    // the executable is not in the current directory
    if(pos != std::string::npos){
        reverse(path.begin(), path.end());
        path_from_exec = path.substr(0, path.length() - pos);
    }
    path_from_exec += POLYEXP_BASE + cli.option('P');

    std::cout << "path: " << path_from_exec << std::endl;
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
    Benchmark_preprocessing b("preprocess", r, false, 0, path_from_exec, timer);
    b.run_benchmark();
}
