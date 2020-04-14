#include "benchmark.hpp"


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

    cli.addOption('r', "100", "number of repetitions");
  
    cliFun.preParse();
    if (!cli.parse()) {return -1;}
    cliFun.postParse();

    int reps = std::stoi(cli.option('r'));

    Test_macro *benchmark = new Test_macro("test_macro", reps, false);
    benchmark->run_benchmark();
}
