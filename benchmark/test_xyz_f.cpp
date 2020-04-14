#include "benchmark.hpp"


class Test_xyz_f : public Benchmark_base_cli {
    public:
        Test_xyz_f(std::string name, int reps, bool convergence, CLIFunctionsVolume &cliFun, bool benchmark_all) : Benchmark_base_cli(name, reps, convergence, cliFun, benchmark_all) {}

        void initialize () {
            box = Polytope_new_box(4,2);
        }
        void reset () {
            // Nothing to reset
        }
        double run () {
            xyz_f(box,0.1,4);
        }
        int get_nr_functions(){
            auto o = dynamic_cast<CLIF_Option<xyz_f_t>*>(cliFun.getOption("xyz_f"));
            return o->fmap.size();
        }
        void select(int s){
            auto o = dynamic_cast<CLIF_Option<xyz_f_t>*>(cliFun.getOption("xyz_f"));
            auto it = o->fmap.begin();
            for (int i = 0; i < s; ++i) {
                it ++;
            }
            selected = it->second;
            name_selected = it->first;
        }
        double run_selected(){
            selected(box,0.1,4);
        }

    private:
        Polytope* box; 
        xyz_f_t selected;
};

int main(int argc, char *argv[]){
    CLI cli(argc,argv,"benchmark");
    CLIFunctionsVolume cliFun(cli);
    
    cli.addOption('r', "100", "number of repetitions");
  
    cliFun.preParse();
    if (!cli.parse()) {return -1;}
    cliFun.postParse();

    int reps = std::stoi(cli.option('r'));

    Test_xyz_f *benchmark = new Test_xyz_f("test_xyz_f", reps, false, cliFun, false);
    benchmark->run_benchmark();
    benchmark = new Test_xyz_f("test_xyz_f", reps, false, cliFun, true);
    benchmark->run_benchmark();
    benchmark = new Test_xyz_f("test_xyz_f", reps, true, cliFun, false);
    benchmark->run_benchmark();
    benchmark = new Test_xyz_f("test_xyz_f", reps, true, cliFun, true);
    benchmark->run_benchmark();
}
