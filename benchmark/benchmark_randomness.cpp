#include "benchmark.hpp"
#include "../src/volume/volume_helper.hpp"

extern "C" {
#include "../src/random/prng.h"
}

enum Value_type {
    random_int,
    random_int_in_range,
    random_double,
    random_double_0_1,
    random_double_normal,
    random_double_in_range
};


class Benchmark_randomness : public Benchmark_base {
public:
    Benchmark_randomness(std::string name,
                            int reps,
                            bool convergence,
                            int warmup_reps,
                            int nr_ints_,
                            int nr_doubles_,
                            Value_type type_) :
        Benchmark_base(name, reps, convergence, warmup_reps, 0.95,0.95), nr_ints(nr_ints_), nr_doubles(nr_doubles_), type(type_)
    {}

protected:
    void initialize () {
        ints = (int*) malloc(nr_ints*sizeof(int));
        doubles = (FT*) malloc(nr_doubles*sizeof(double));
        prng_init();
    }
    
    void reset () {}

    double run () {
        switch(type) {
            case random_int:
                for(int i = 0; i < nr_ints; ++i) {
                    ints[i] = prng_get_random_int();
                };
                break;
            case random_int_in_range:
                for(int i = 0; i < nr_ints; ++i) {
                    ints[i] = prng_get_random_int_in_range(-63, 64);
                };
                break;
            case random_double:
                for(int i = 0; i < nr_doubles; ++i) {
                    doubles[i] = prng_get_random_double();
                };
                break;
            case random_double_0_1:
                for(int i = 0; i < nr_doubles; ++i) {
                    doubles[i] = prng_get_random_double_0_1();
                };
                break;
            case random_double_normal:
                for(int i = 0; i < nr_doubles; ++i) {
                    doubles[i] = prng_get_random_double_normal();
                };
                break;
            case random_double_in_range:
                for(int i = 0; i < nr_doubles; ++i) {
                    doubles[i] = prng_get_random_double_in_range(-63.0, 64.0);
                };
                break;
            default:
                assert(1 == 0 && "Incomplete switch for emum: Type not covered");
        }
        return 0;
    }
    
    void finalize() {
	    Benchmark_base::finalize();
        free(ints);
        free(doubles);
	}
    
private:
    int *ints;
    FT *doubles;
    int nr_ints;
    int nr_doubles;
    Value_type type;
};

int main(int argc, char *argv[]){
    CLI cli(argc,argv,"benchmark");
    CLIFunctionsVolume cliFun(cli);

    int r = 100;
    cliFun.claimOpt('b',"Benchmarking configuration");
    cliFun.add(new CLIF_OptionNumber<int>(&r,'b',"r","100", 1, 1000000)); // 1 page of ints

    int nr_ints = 2048;
    cliFun.add(new CLIF_OptionNumber<int>(&nr_ints,'b',"i","100", 1, 1048576)); // 1 page of ints

    int nr_doubles = 1024;
    cliFun.add(new CLIF_OptionNumber<int>(&nr_doubles,'b',"d","100", 1, 524288)); // 1 page of doubles

    Value_type v_t = random_int;
    cliFun.add(new CLIF_Option<Value_type>(&v_t,'b',"rand_val_t","random_int", {
                        {"random_int",{random_int, "prng_get_random_int"}},
					    {"random_int_in_range",{random_int_in_range, "prng_get_random_int_in_range"}},
					    {"random_double",{random_double, "prng_get_random_double"}},
					    {"random_double_0_1",{random_double_0_1, "prng_get_random_double_0_1"}},
					    {"random_double_normal",{random_double_normal, "prng_get_random_double_normal"}},
					    {"random_double_in_range",{random_double_in_range, "prng_get_random_double_in_range"}} }));
    
    cliFun.preParse();
    if (!cli.parse()) {return -1;}
    cliFun.postParse();
    
    Benchmark_randomness b("randomness", r, true, 10, nr_ints, nr_doubles, v_t);
    b.run_benchmark();
}
