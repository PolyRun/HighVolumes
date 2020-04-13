#include <iostream>
#include <cmath>
#include <limits>
#include <algorithm>
#include <vector>
#include <string>

#include "../src/util/timer.hpp"
#include "../src/volume/volume_helper.hpp"
#include "../polyvest/vol.h"

#include "../src/util/cli.hpp"
#include "../src/util/cli_functions.hpp"

extern "C" { // must be included C stlye
#include "../src/random/prng.h"
}

#include "benchmark.hpp"

#define MACRO_REPS 100
#define MINI_REPS 100

int macro_benchmarks_count = 0;
std::vector<Benchmark_base*> macro_benchmarks;

void add_macro_function(Benchmark_base &b);


int mini_benchmarks_count = 0;
std::vector<Benchmark_base_cli*> mini_benchmarks;

void add_mini_function(Benchmark_base_cli &b);

void macro_benchmark_polyvest(Timer &timer, double &min_time, double &max_time, double &mean_time, double& std_dev) {
    /* Prepare input */
    Polytope *P;
    read_polyvest_p("../polyvest/examples/cube_2", &P);
    
    int n = P->n;
    int m = P->m;
    
    vol::Polyvest_p Q(m, n);
  
    /*vol::Polyvest_p Q(m, n);
    polyvest_convert(P, &Q);*/
    double total_time = 0;
    double measured_times[MACRO_REPS];
    for (int i = 0; i < MACRO_REPS; ++i) {

        // Polyvest changes polytope while processing it
        polyvest_convert(P, &Q);

        // Run benchmark
        timer.start();
        Q.Preprocess();
        Q.EstimateVol(1600);
        timer.stop();
        
        measured_times[i] = timer.millisecs();            
        total_time += measured_times[i];
        min_time = std::min(min_time, measured_times[i]);
        max_time = std::max(max_time, measured_times[i]);

    }
    
    mean_time = total_time/MACRO_REPS;
    for (int i = 0; i < MACRO_REPS; ++i) {
        std_dev += pow(measured_times[i] - mean_time, 2.0);
    }
    std_dev = sqrt(std_dev/MACRO_REPS);
}

/**
 * \brief Benchmarks for smaller functions in order to compare different versions
 **/
void run_mini_benchmarks(){
    Timer timer;

    for(int i = 0; i < mini_benchmarks_count; ++i) {
        Benchmark_base_cli *b = mini_benchmarks[i];
        if (b->get_benchmark_all()) {
            for(int j = 0; j < b->get_nr_functions(); j++) {
                double min_time = std::numeric_limits<double>::max();
                double max_time = -1;
                double mean_time;
                double std_dev = 0.0;
                double total_time = 0;
                double measured_times[MINI_REPS];

                // Initialize
                b->select(j);
                b->initialize();

                for (int i = 0; i < MINI_REPS; ++i) {

                    // Run benchmark
                    timer.start();
                    b->run_selected();
                    timer.stop();
                    
                    measured_times[i] = timer.millisecs();            
                    total_time += measured_times[i];
                    min_time = std::min(min_time, measured_times[i]);
                    max_time = std::max(max_time, measured_times[i]);

                    // Reset
                    b->reset();

                }
                
                mean_time = total_time/MINI_REPS;
                for (int i = 0; i < MINI_REPS; ++i) {
                    std_dev += pow(measured_times[i] - mean_time, 2.0);
                }
                std_dev = sqrt(std_dev/MINI_REPS);
                
                std::cout << "name: "<< b->get_name_selected() << ", mean: " << mean_time << ", min: " << min_time << ", max: " << max_time << ", std dev: " << std_dev << std::endl;
            }
        } else {
            double min_time = std::numeric_limits<double>::max();
            double max_time = -1;
            double mean_time;
            double std_dev = 0.0;
            double total_time = 0;
            double measured_times[MINI_REPS];

            // Initialize
            b->initialize();

            for (int i = 0; i < MINI_REPS; ++i) {

                // Run benchmark
                timer.start();
                b->run();
                timer.stop();
                
                measured_times[i] = timer.millisecs();            
                total_time += measured_times[i];
                min_time = std::min(min_time, measured_times[i]);
                max_time = std::max(max_time, measured_times[i]);

                // Reset
                b->reset();

            }
            
            mean_time = total_time/MINI_REPS;
            for (int i = 0; i < MINI_REPS; ++i) {
                std_dev += pow(measured_times[i] - mean_time, 2.0);
            }
            std_dev = sqrt(std_dev/MINI_REPS);
            
            std::cout << "name: "<< b->get_name() << ", mean: " << mean_time << ", min: " << min_time << ", max: " << max_time << ", std dev: " << std_dev << std::endl;
        }
        
    }
    
}

/**
 * \brief Benchmarks for large functions in order to generate plot data
 **/
void run_macro_benchmarks(){
    Timer timer;

    for(int i = 0; i < macro_benchmarks_count; ++i) {
        Benchmark_base *b = macro_benchmarks[i];
        double min_time = std::numeric_limits<double>::max();
        double max_time = -1;
        double mean_time;
        double std_dev = 0.0;
        double total_time = 0;
        double measured_times[MACRO_REPS];

        // Initialize
        b->initialize();

        for (int i = 0; i < MACRO_REPS; ++i) {

            // Run benchmark
            timer.start();
            b->run();
            timer.stop();
            
            measured_times[i] = timer.millisecs();            
            total_time += measured_times[i];
            min_time = std::min(min_time, measured_times[i]);
            max_time = std::max(max_time, measured_times[i]);

            // Reset
            b->reset();

        }
        
        mean_time = total_time/MACRO_REPS;
        for (int i = 0; i < MACRO_REPS; ++i) {
            std_dev += pow(measured_times[i] - mean_time, 2.0);
        }
        std_dev = sqrt(std_dev/MACRO_REPS);
        
        std::cout << "name: "<< b->get_name() << ", mean: " << mean_time << ", min: " << min_time << ", max: " << max_time << ", std dev: " << std_dev << std::endl;
    }

}

void add_macro_function(Benchmark_base &b) {
    macro_benchmarks.push_back(&b);
    macro_benchmarks_count++;
}

void add_mini_function(Benchmark_base_cli &b) {
    mini_benchmarks.push_back(&b);
    mini_benchmarks_count++;
}


void add_macro_benchmark_functions(){
    add_macro_function(*(new Macro_benchmark_test("macro_benchmark_test_name")));
    //add_macro_function(macro_benchmark_polyvest, "macro_benchmark_polyvest");
}
    
void add_mini_benchmark_functions(CLIFunctionsVolume &cliFun){
    add_mini_function(*(new Mini_benchmark_xyz_f("mini_benchmark_xyz_f", cliFun, true)));
    add_mini_function(*(new Mini_benchmark_xyz_f("mini_benchmark_xyz_f", cliFun, false)));
}


int main(int argc, char *argv[]){
    CLI cli(argc,argv,"benchmark");
    CLIFunctionsVolume cliFun(cli);
  
    cliFun.preParse();
    if (!cli.parse()) {return -1;}
    cliFun.postParse();


    std::cout << "===== Starting benchmarks =====\n" << std::endl;

    /* Add macro-benchmark functions */
    add_macro_benchmark_functions();
    
    /* Add mini-benchmark functions */
    add_mini_benchmark_functions(cliFun);


    /* Run function-benchmarks */
    std::cout << "\n== Run mini-benchmarks with " << MINI_REPS << " repetitions each ==\n" << std::endl;
    run_mini_benchmarks();

    /* Run macro-benchmarks */
    std::cout << "\n== Run macro-benchmarks with " << MACRO_REPS << " repetitions each ==\n" << std::endl;
    run_macro_benchmarks();

    std::cout << "\n===== Finished benchmarks =====" << std::endl;
}
