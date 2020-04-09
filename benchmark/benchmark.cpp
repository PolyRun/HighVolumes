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


typedef void(*macro_fun)(Timer &timer, double &min_time, double &max_time, double &mean_time, double &std_dev);

typedef void(*mini_fun)(Timer &selected_timer, double &selected_min_time, double &selected_max_time, double &selected_mean_time, double &selected_std_dev, bool benchmark_all, CLIFunctionsVolume &cliFun);

int macro_functions_count = 0;
std::vector<macro_fun> macro_functions;
std::vector<std::string> macro_functions_names;

void add_macro_function(macro_fun f, std::string name);


int mini_functions_count = 0;
std::vector<mini_fun> mini_functions;
std::vector<std::string> mini_functions_names;
std::vector<bool> mini_functions_benchmark_all;

void add_mini_function(mini_fun f, std::string name);

void mini_benchmark_xyz_f(Timer &timer, double &selected_min_time, double &selected_max_time, double &selected_mean_time, double &selected_std_dev, bool benchmark_all, CLIFunctionsVolume &cliFun) {
    /* Prepare input */
    // Generate new polytope box, 4 dim, 2 radius
    Polytope* box = Polytope_new_box(4,2);

    /* Benchmark selected version */
    double selected_total_time = 0;
    double selected_measured_times[MINI_REPS];
    for (int i = 0; i < MINI_REPS; ++i) {

        // Run benchmark
        timer.start();
        xyz_f(box,0.1,4);
        timer.stop();
        
        selected_measured_times[i] = timer.millisecs();            
        selected_total_time += selected_measured_times[i];
        selected_min_time = std::min(selected_min_time, selected_measured_times[i]);
        selected_max_time = std::max(selected_max_time, selected_measured_times[i]);

    }
    
    selected_mean_time = selected_total_time/MINI_REPS;
    for (int i = 0; i < MINI_REPS; ++i) {
        selected_std_dev += pow(selected_measured_times[i] - selected_mean_time, 2.0);
    }
    selected_std_dev = sqrt(selected_std_dev/MINI_REPS);

    /* Benchmark all versions */
    if(benchmark_all) {
        auto o = dynamic_cast<CLIF_Option<xyz_f_t>*>(cliFun.getOption("xyz_f"));
        for(auto it : o->fmap) {
            double opt_min_time = std::numeric_limits<double>::max();
            double opt_max_time = -1;
            double opt_mean_time;
            double opt_std_dev = 0.0;

            double opt_total_time = 0;
            double opt_measured_times[MINI_REPS];
            for (int i = 0; i < MINI_REPS; ++i) {

                // Run benchmark
                timer.start();
                it.second(box,0.1,4);
                timer.stop();
                
                opt_measured_times[i] = timer.millisecs();            
                opt_total_time += opt_measured_times[i];
                opt_min_time = std::min(opt_min_time, opt_measured_times[i]);
                opt_max_time = std::max(opt_max_time, opt_measured_times[i]);

            }
            
            opt_mean_time = opt_total_time/MINI_REPS;
            for (int i = 0; i < MINI_REPS; ++i) {
                opt_std_dev += pow(opt_measured_times[i] - opt_mean_time, 2.0);
            }
            opt_std_dev = sqrt(opt_std_dev/MINI_REPS);

            std::cout << "name(opt): "<<it.first << ", mean: " << opt_mean_time << ", min: " << opt_min_time << ", max: " << opt_max_time << ", std dev: " << opt_std_dev << std::endl;
        }
    }
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
 * \brief Benchmarks for smaller functions in order to compare different versions
 **/
void run_mini_benchmarks(CLIFunctionsVolume &cliFun){
    Timer timer;

    for(int i = 0; i < mini_functions_count; ++i) {
        mini_fun f = mini_functions[i];
        double selected_min_time = std::numeric_limits<double>::max();
        double selected_max_time = -1;
        double selected_mean_time;
        double selected_std_dev = 0.0;

        // Benchmark function
        f(timer, selected_min_time, selected_max_time, selected_mean_time, selected_std_dev, mini_functions_benchmark_all[i], cliFun);
        
        std::cout << "name(sel): "<<mini_functions_names[i] << ", mean: " << selected_mean_time << ", min: " << selected_min_time << ", max: " << selected_max_time << ", std dev: " << selected_std_dev << std::endl;
    }
    
}

/**
 * \brief Benchmarks for large functions in order to generate plot data
 **/
void run_macro_benchmarks(){
    Timer timer;

    for(int i = 0; i < macro_functions_count; ++i) {
        macro_fun f = macro_functions[i];
        double min_time = std::numeric_limits<double>::max();
        double max_time = -1;
        double mean_time;
        double std_dev = 0.0;

        // Benchmark function
        f(timer, min_time, max_time, mean_time, std_dev);
        
        std::cout << "name: "<<macro_functions_names[i] << ", mean: " << mean_time << ", min: " << min_time << ", max: " << max_time << ", std dev: " << std_dev << std::endl;
    }

}

void add_macro_function(macro_fun f, std::string name) {
    macro_functions.push_back(f);
    macro_functions_names.push_back(name);
    macro_functions_count++;
}

void add_mini_function(mini_fun f, std::string name, bool benchmark_all) {
    mini_functions.push_back(f);
    mini_functions_names.push_back(name);
    mini_functions_benchmark_all.push_back(benchmark_all);
    mini_functions_count++;
}


void add_macro_benchmark_functions(){
    add_macro_function(macro_benchmark_test, "macro_benchmark_test");
    add_macro_function(macro_benchmark_polyvest, "macro_benchmark_polyvest");    
}
    
void add_mini_benchmark_functions(){
    add_mini_function(mini_benchmark_xyz_f, "mini_benchmark_xyz_f", true);
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
    add_mini_benchmark_functions();


    /* Run function-benchmarks */
    std::cout << "\n== Run mini-benchmarks with " << MINI_REPS << " repetitions each ==\n" << std::endl;
    run_mini_benchmarks(cliFun);

    /* Run macro-benchmarks */
    std::cout << "\n== Run macro-benchmarks with " << MACRO_REPS << " repetitions each ==\n" << std::endl;
    run_macro_benchmarks();

    std::cout << "\n===== Finished benchmarks =====" << std::endl;
}
