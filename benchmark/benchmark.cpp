#include <iostream>
#include <cmath>
#include <limits>
#include <algorithm>
#include <vector>
#include <string>

#include "../src/util/timer.hpp"
#include "../src/poly/volume_helper.hpp"
#include "../polyvest/vol.h"

extern "C" { // must be included C stlye
#include "../src/random/prng.h"
}

//#include "benchmark.hpp"

#define MACRO_REPS 100
#define MINI_REPS 100

int read_polyvest_p(string filename, Polytope **P){

    ifstream file;
    file.open(filename);


    if (!file.is_open()){
        printf("failed to read polytope");
        return 1;
    }

    int n, m;
    file >> m >> n;

    *P = Polytope_new(n, m);

    FT num;
    for (int i = 0; i < m; i++){
        file >> num;
        Polytope_set_b(*P, i, num);
        for (int j = 0; j < n; j++){
            file >> num;
            Polytope_set_a(*P, i, j, num);
        }
    }

    return 0;
}

void polyvest_convert(Polytope *P, vol::Polyvest_p *Q){

    int n = P->n;
    int m = P->m;  

    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            Q->matA(Polytope_get_a(P, i, j), i, j);
        }
        Q->vecb(Polytope_get_b(P, i), i);
    }
  
}


typedef void(*fun)(Timer &timer, double &min_time, double &max_time, double &mean_time, double &std_dev);

int macro_functions_count = 0;
std::vector<fun> macro_functions;
std::vector<std::string> macro_functions_names;

void add_macro_function(fun f, std::string name);


/**
 * \brief Benchmarks for smaller functions in order to compare different versions
 **/
void run_mini_benchmarks(){
    // TODO
}

void macro_benchmark_test(Timer &timer, double &min_time, double &max_time, double &mean_time, double& std_dev) {
    /* Prepare input */
    
    double total_time = 0;
    double measured_times[MACRO_REPS];
    for (int i = 0; i < MACRO_REPS; ++i) {

        // Run benchmark
        timer.start();
        prng_get_random_int_in_range(1, 2020);
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
 * \brief Benchmarks for large functions in order to generate plot data
 **/
void run_macro_benchmarks(){
    Timer timer;

    // Prepare input

    for(int i = 0; i < macro_functions_count; ++i) {
        fun f = macro_functions[i];
        double min_time = std::numeric_limits<double>::max();
        double max_time = -1;
        double mean_time;
        double std_dev = 0.0;

        // Benchmark function
        f(timer, min_time, max_time, mean_time, std_dev);
        
        std::cout << "name: "<<macro_functions_names[i] << ", mean: " << mean_time << ", min: " << min_time << ", max: " << max_time << ", std dev: " << std_dev << std::endl;
    }
    
}

void add_macro_function(fun f, std::string name) {
    macro_functions.push_back(f);
    macro_functions_names.push_back(name);
    macro_functions_count++;
}


void add_macro_benchmark_functions(){
    add_macro_function(macro_benchmark_test, "macro_benchmark_test");
    add_macro_function(macro_benchmark_polyvest, "macro_benchmark_polyvest");    
}
    
void add_mini_benchmark_functions(){
    // TODO
}


int main(int argc, char *argv[]){
    std::cout << "===== Starting benchmarks =====\n" << std::endl;

    /* Add macro-benchmark functions */
    add_macro_benchmark_functions();
    
    /* Add mini-benchmark functions */
    add_mini_benchmark_functions();


    /* Run function-benchmarks */
    run_mini_benchmarks();
    std::cout << "\n== Run mini-benchmarks with " << MINI_REPS << " repetitions each ==\n" << std::endl;

    /* Run macro-benchmarks */
    std::cout << "\n== Run macro-benchmarks with " << MACRO_REPS << " repetitions each ==\n" << std::endl;
    run_macro_benchmarks();

    std::cout << "\n===== Finished benchmarks =====" << std::endl;
}
