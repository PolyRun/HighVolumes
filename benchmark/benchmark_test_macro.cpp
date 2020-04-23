#include "benchmark.hpp"


class Test_macro : public Benchmark_base {
    public:
        Test_macro(std::string name, int reps, bool convergence, int warmup_reps) : Benchmark_base(name, reps, convergence, warmup_reps) {}

    protected:
        void initialize () {
            std::cout << "initializing macro benchmark test" << std::endl;
        }
        void reset () {
            std::cout << "resetting macro benchmark test" << std::endl;
        }
        double run () {
            std::cout << "running macro benchmark test" << std::endl;
        }
};

int main(int argc, char *argv[]){
    CLI cli(argc,argv,"benchmark");
    CLIFunctions cliFun(cli);

    int r = 100;
    cliFun.claimOpt('b',"Benchmarking configuration");
    cliFun.add(new CLIF_OptionNumber<int>(&r,'b',"r","100", 1, 100000));
  
    cliFun.preParse();
    if (!cli.parse()) {return -1;}
    cliFun.postParse();

    Test_macro *benchmark = new Test_macro("test_macro", r, false, 0);
    benchmark->run_benchmark();
}
