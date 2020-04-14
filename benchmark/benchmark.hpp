#include <iostream>
#include <string>
#include <cmath>
#include <limits>
#include <algorithm>

#include "../src/util/cli.hpp"
#include "../src/util/cli_functions.hpp"
#include "../src/util/timer.hpp"

#ifndef BENCHMARK_H
#define BENCHMARK_H


/* Base classes */

/**
 * Base class for benchmarking functions (abstract)
 **/
class Benchmark_base {
    public:
        Benchmark_base(std::string name_, int reps_, bool convergence_) : name(name_), reps(reps_), convergence(convergence_){}

        /**
         * \brief Actuall performs the benchmark
         **/
        virtual double run_benchmark() {
            double min_time = std::numeric_limits<double>::max();
            double max_time = -1;
            double mean_time;
            double std_dev = 0.0;
            double total_time = 0;
            double measured_times[reps];

            // Initialize
            initialize();

            for (int i = 0; i < reps; ++i) {

                // Run benchmark
                timer.start();
                double result = run();
                timer.stop();

                if (convergence) {
                    if (i == 0) {
                        std::cout << "Convergence Step(" << i << ") - Result: " << result << std::endl;
                    } else {
                        std::cout << "Convergence Step(" << i << ") - Result: " << result << ", Diff: " << result-last_result << std::endl;
                    }
                    last_result = result;
                }
                
                measured_times[i] = timer.millisecs();            
                total_time += measured_times[i];
                min_time = std::min(min_time, measured_times[i]);
                max_time = std::max(max_time, measured_times[i]);

                // Reset
                reset();

            }
            
            mean_time = total_time/reps;
            for (int i = 0; i < reps; ++i) {
                std_dev += pow(measured_times[i] - mean_time, 2.0);
            }
            std_dev = sqrt(std_dev/reps);
            
            std::cout << "name: "<< name << ", mean: " << mean_time << ", min: " << min_time << ", max: " << max_time << ", std dev: " << std_dev << std::endl;
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

        std::string name; // Name of the benchmark that is displayed in output
        int reps; // Number of repetitions in benchmark
        bool convergence; // Optional convergence output
        double last_result; // Value in last step;
        Timer timer;
};

/**
 * Base class for benchmarking functions that use cli (abstract)
 **/
class Benchmark_base_cli : public Benchmark_base{
    public:
        Benchmark_base_cli(std::string name_, int reps_, bool convergence_, CLIFunctions &cliFun_, bool benchmark_all_) : Benchmark_base(name_, reps_, convergence_), cliFun(cliFun_), benchmark_all(benchmark_all_){
        }

        virtual double run_benchmark() {
            for (int j = 0; j < get_nr_functions(); ++j){
                double min_time = std::numeric_limits<double>::max();
                double max_time = -1;
                double mean_time;
                double std_dev = 0.0;
                double total_time = 0;
                double measured_times[reps];

                // Initialize
                initialize();
                select(j);

                for (int i = 0; i < reps; ++i) {
                    double result;
                    // Run benchmark
                    if (benchmark_all) {
                        timer.start();
                        result = run_selected();
                        timer.stop();
                    } else {
                        timer.start();
                        result = run();
                        timer.stop();
                    }

                    if (convergence) {
                        if (i == 0) {
                            std::cout << "Convergence Step(" << i << ") - Result: " << result << std::endl;
                        } else {
                            std::cout << "Convergence Step(" << i << ") - Result: " << result << ", Diff: " << result-last_result << std::endl;
                        }
                        last_result = result;
                    }
                    
                    measured_times[i] = timer.millisecs();            
                    total_time += measured_times[i];
                    min_time = std::min(min_time, measured_times[i]);
                    max_time = std::max(max_time, measured_times[i]);

                    // Reset
                    reset();

                }
                
                mean_time = total_time/reps;
                for (int i = 0; i < reps; ++i) {
                    std_dev += pow(measured_times[i] - mean_time, 2.0);
                }
                std_dev = sqrt(std_dev/reps);
                
                if (benchmark_all) {
                    std::cout << "name: "<< name_selected << ", mean: " << mean_time << ", min: " << min_time << ", max: " << max_time << ", std dev: " << std_dev << std::endl;
                } else {
                    std::cout << "name: "<< name << ", mean: " << mean_time << ", min: " << min_time << ", max: " << max_time << ", std dev: " << std_dev << std::endl;
                    break;
                }
            }
        }

    protected:
        /**
         * Returns the number of functions that are available
         **/
        virtual int get_nr_functions() = 0;

        /**
         * Selects the function that should be run
         **/
        virtual void select(int s) = 0;

        /**
         * Runs the selected function (only function call to prevent overhead)
         **/
        virtual double run_selected() = 0;

        CLIFunctions cliFun;
        bool benchmark_all; // If set, all available functions will be benchmarked, otherwise only the one selected by cli
        std::string name_selected; // Name of the currently selected function
};

#endif // BENCHMARK_H
