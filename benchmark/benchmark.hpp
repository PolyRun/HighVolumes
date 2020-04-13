
#include <iostream>
#include <string>

/* Base class includes */
#include "../src/util/cli.hpp"
#include "../src/util/cli_functions.hpp"

/* Additional classes includes */
#include "../src/volume/volume_helper.hpp"
#include "../polyvest/vol.h"

#ifndef BENCHMARK_H
#define BENCHMARK_H


/* Base classes */

/**
 * Base class for benchmarking functions (abstract)
 **/
class Benchmark_base {
    public:
        Benchmark_base(std::string name_) : name(name_){}

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

        /**
         * Setter for benchmark name
         **/
        void set_name(std::string name_) {
            name = name_;
        }

        /**
         * Getter for benchmark name
         **/
        std::string get_name(){
            return name;
        }

    protected:
        std::string name; // Name of the benchmark that is displayed in output
};

/**
 * Base class for benchmarking functions that use cli (abstract)
 **/
class Benchmark_base_cli : public Benchmark_base{
    public:
        Benchmark_base_cli(std::string name_, CLIFunctionsVolume &cliFun_, bool benchmark_all_) : Benchmark_base(name_), cliFun(cliFun_), benchmark_all(benchmark_all_){
        }

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

        /**
         * Setter for benchmark_all
         **/
        void set_benchmark_all(bool b) {
            benchmark_all = b;
        }

        /**
         * Getter for benchmark_all
         **/
        bool get_benchmark_all(){
            return benchmark_all;
        }
        
        /**
         * Getter for the name of the currently selected function
         **/
        std::string get_name_selected(){
            return name_selected;
        }

    protected:
        CLIFunctionsVolume cliFun;
        bool benchmark_all; // If set, all available functions will be benchmarked, otherwise only the one selected by cli
        std::string name_selected; // Name of the currently selected function
};


/* Macro benchmarks */

class Macro_benchmark_test : public Benchmark_base {
    public:
        Macro_benchmark_test(std::string name) : Benchmark_base(name) {}

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

/* Mini benchmarks */

class Mini_benchmark_xyz_f : public Benchmark_base_cli {
    public:
        Mini_benchmark_xyz_f(std::string name, CLIFunctionsVolume &cliFun, bool benchmark_all) : Benchmark_base_cli(name, cliFun, benchmark_all) {}

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


#endif // BENCHMARK_H
