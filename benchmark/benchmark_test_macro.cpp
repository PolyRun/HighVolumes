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

    cli.addOption('r', "100", "number of repetitions");
  
    cliFun.preParse();
    if (!cli.parse()) {return -1;}
    cliFun.postParse();

    int reps = std::stoi(cli.option('r'));

    Test_macro *benchmark = new Test_macro("test_macro", reps, false, 0);
    benchmark->run_benchmark();
}
