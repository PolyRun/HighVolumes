#include <iostream>
#include <string>
#include <cmath>
#include <limits>
#include <algorithm>
#include <vector>
#include <omp.h>

#include "../src/util/cli.hpp"
#include "../src/util/cli_functions.hpp"
#include "../src/util/timer.hpp"
#include "../src/util/tsc_x86.h"

#ifndef BENCHMARK_H
#define BENCHMARK_H


/* Base classes */

/**
 * Base class for benchmarking functions (abstract)
 **/
class Benchmark_base {
    public:
    Benchmark_base(std::string name_, int reps_, bool convergence_, int warmup_reps_, double time_ci_alpha_, double results_ci_alpha_)
	    : name(name_), reps(reps_), convergence(convergence_), warmup_reps(warmup_reps_), time_ci_alpha(time_ci_alpha_), results_ci_alpha(results_ci_alpha_) {
        timer = new Tsc();
    }

    // overloaded constructer taking timer argument
    Benchmark_base(std::string name_, int reps_, bool convergence_, int warmup_reps_, double time_ci_alpha_, double results_ci_alpha_, Timer_generic *t)
	    : name(name_), reps(reps_), convergence(convergence_), warmup_reps(warmup_reps_), time_ci_alpha(time_ci_alpha_), results_ci_alpha(results_ci_alpha_), timer(t) {
        // empty
    }

        /**
         * \brief Actuall performs the benchmark
         **/
        virtual void run_benchmark() {
            run_benchmark_go();
	    //#pragma omp parallel
	    //{
	    //   int tid = omp_get_thread_num();
	    //   #pragma omp barrier
	    //   if(tid == 0) {
	    //      int nthreads = omp_get_num_threads();
	    //      printf("There are %d threads, running on tid %d.\n",nthreads,tid);
	    //      run_benchmark_go();
	    //   }
	    //}
	}
	void run_benchmark_go() {
	    double min_time = std::numeric_limits<double>::max();
            double max_time = -1;
            double mean_time;
            double std_dev = 0.0;
            double total_time = 0;
	    std::vector<double> measured_times(reps);
            
	    std::vector<double> results(reps);
            double results_sum = 0.0;

            // Initialize
            initialize();

	    std::cout << "Benchmark Warmup ("<< warmup_reps <<")\n";
            for (int i = 0; i < warmup_reps; ++i) {
                run();
		reset();
            }

	    std::cout << "Benchmark Measurement ("<< reps <<")\n";
            for (int i = 0; i < reps; ++i) {

                // Run benchmark
                timer->start();
                double result = run();
                timer->stop();

                results[i] = result;
                results_sum += result;

                measured_times[i] = timer->get_time();            
                total_time += measured_times[i];
                min_time = std::min(min_time, measured_times[i]);
                max_time = std::max(max_time, measured_times[i]);

                // Reset
                reset();

            }
            
	    // process time:
            mean_time = total_time/reps;
            for (int i = 0; i < reps; ++i) {
                std_dev += pow(measured_times[i] - mean_time, 2.0);
            }
            std_dev = sqrt(std_dev/reps);

	    // Run performance_counter and free memory
	    finalize();

	    // Calculate confidence intervals:
	    double time_ci_low, time_ci_high, time_ci_median;
	    {
	       std::sort(measured_times.begin(), measured_times.end());
	       size_t s = measured_times.size();
               double time_ci_a1 = (1.0-0.95)/2.0;
	       double time_ci_a2 = time_ci_a1 + time_ci_alpha;
	       time_ci_low    = measured_times[s*time_ci_a1];
	       time_ci_high   = measured_times[s*time_ci_a2];
	       time_ci_median = measured_times[s*0.5];
	    }
	
	    if(pc_flops > 0 || pc_bytes > 0 ) {
	       std::cout << "Avg flops/cycle: " << pc_flops/mean_time << " (CI: " << pc_flops/time_ci_low << " - " << pc_flops/time_ci_high << ")\n";
	       std::cout << "Avg bytes/cycle: " << pc_bytes/mean_time << " (CI: " << pc_bytes/time_ci_low << " - " << pc_bytes/time_ci_high << ")\n";
	    }
            
            

            // Dictionary output
            std::cout << "{'time': {";
            std::cout << "'name_t': '"<< name << "', 'mean': '" << mean_time << "', 'min': '" << min_time << "', 'max': '" << max_time << "', 'std_dev': '" << std_dev << "'";
	    
            std::cout << ", 'ci_low': '"<< time_ci_low <<"'";
	    std::cout << ", 'ci_high': '"<< time_ci_high <<"'";
	    std::cout << ", 'ci_median': '"<< time_ci_median <<"'";
	    std::cout << ", 'ci_alpha': '"<< time_ci_alpha <<"'";

            std::cout << "}, 'convergence': {";
	        // process results:
            if(convergence) {
                double results_mean = results_sum/reps;
                double results_min = results[0];
                double results_max = results[0];
                double results_std_dev = 0;
    
                for (int i = 0; i < reps; ++i) {
                    results_std_dev += pow(results[i] - results_mean, 2.0);
                    results_max = std::max(results_max, results[i]);
                    results_min = std::min(results_min, results[i]);
                }
                results_std_dev = sqrt(results_std_dev/reps);
                    
                std::cout << "'name_c': '"<< name << "', 'mean': '" << results_mean << "', 'min': '" << results_min << "', 'max': '" << results_max << "', 'std_dev': '" << results_std_dev << "'";
            }

	    // Calculate confidence intervals:
	    double results_ci_low, results_ci_high, results_ci_median;
	    {
	       std::sort(results.begin(), results.end());
	       size_t s = results.size();
               double results_ci_a1 = (1.0-0.95)/2.0;
	       double results_ci_a2 = results_ci_a1 + results_ci_alpha;
	       results_ci_low    = results[s*results_ci_a1];
	       results_ci_high   = results[s*results_ci_a2];
	       results_ci_median = results[s*0.5];
	    }
	    std::cout << ", 'ci_low': '"<< results_ci_low <<"'";
	    std::cout << ", 'ci_high': '"<< results_ci_high <<"'";
	    std::cout << ", 'ci_median': '"<< results_ci_median <<"'";
	    std::cout << ", 'ci_alpha': '"<< results_ci_alpha <<"'";

	    std::cout << "}, 'performance_counter': {";
	    std::cout << "'flops': " << pc_flops << ", 'bytes': " << pc_bytes;
            std::cout << "}}" << std::endl;
        }

    protected:

        /**
         * Initializes all data that is needed in order to run the function (e.g. input data)
         **/
        virtual void initialize() = 0;

        /**
         * Prepares everything for the function to be rerun.
         **/
        virtual void reset() = 0;

        /**
         * Runs the function (only function call to prevent overhead)
         **/
        virtual double run() = 0;
        
	/*
	 * Called at end, just before measurements are printed
	 * Use to free up memory, and run performance counter (flops / bytes)
	 **/
	virtual void finalize() {
	    std::cout << "Finalize (defalut).\n";
	    pc_flops = 0;
	    pc_bytes = 0;
	}

        std::string name; // Name of the benchmark that is displayed in output
        int reps; // Number of repetitions in benchmark
        bool convergence; // Optional convergence output
        double last_result; // Value in last step;
        Timer_generic *timer;
        int warmup_reps;

	// performance_counter results:
	size_t pc_flops = 0;
	size_t pc_bytes = 0;

        // alphas for confidence intervals
	double time_ci_alpha = 0.95;
	double results_ci_alpha = 0.95;
};



#endif // BENCHMARK_H
