#include "benchmark.hpp"

#define REPS 100

class Test_macro : public Benchmark_base {
    public:
        Test_macro(std::string name, int reps, bool convergence) : Benchmark_base(name, reps, convergence) {}

        void initialize () {
            std::cout << "initializing macro benchmark test" << endl;
        }
        void reset () {
            std::cout << "resetting macro benchmark test" << endl;
        }
        double run () {
            std::cout << "running macro benchmark test" << endl;
        }
};

int main(int argc, char *argv[]){
    CLI cli(argc,argv,"benchmark");
    CLIFunctionsVolume cliFun(cli);
  
    cliFun.preParse();
    if (!cli.parse()) {return -1;}
    cliFun.postParse();

    Test_macro *benchmark = new Test_macro("test_macro", REPS, false);
    benchmark->run_benchmark();
}
