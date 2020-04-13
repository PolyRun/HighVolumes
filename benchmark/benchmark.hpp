#include <string>
#include <iostream>

#ifndef BENCHMARK_H
#define BENCHMARK_H

class Benchmark_base {
    public:
        Benchmark_base(std::string name_) : name(name_){}

        virtual void initialize() = 0;
        virtual void reset() = 0;
        virtual void run() = 0;

        void set_name(std::string name_) {
            name = name_;
        }
        std::string get_name(){
            return name;
        }

    protected:
        std::string name;
};

class Benchmark_base_cli : public Benchmark_base{
    public:
        Benchmark_base_cli(std::string name_, CLIFunctionsVolume &cliFun_, bool benchmark_all_) : Benchmark_base(name_), cliFun(cliFun_), benchmark_all(benchmark_all_){
        }

        virtual int get_nr_functions() = 0;
        virtual void select(int s) = 0;
        virtual void run_selected() = 0;

        void set_benchmark_all(bool b) {
            benchmark_all = b;
        }
        bool get_benchmark_all(){
            return benchmark_all;
        }
        
        std::string get_name_selected(){
            return name_selected;
        }

    protected:
        bool benchmark_all;
        CLIFunctionsVolume cliFun;
        std::string name_selected;
};

class Macro_benchmark_test : public Benchmark_base {
    public:
        Macro_benchmark_test(std::string name) : Benchmark_base(name) {}

        void initialize () {
            std::cout << "initializing macro benchmark test" << endl;
        }
        void reset () {
            std::cout << "resetting macro benchmark test" << endl;
        }
        void run () {
            std::cout << "running macro benchmark test" << endl;
        }
};

class Mini_benchmark_xyz_f : public Benchmark_base_cli {
    public:
        Mini_benchmark_xyz_f(std::string name, CLIFunctionsVolume &cliFun, bool benchmark_all) : Benchmark_base_cli(name, cliFun, benchmark_all) {}

        void initialize () {
            box = Polytope_new_box(4,2);
        }
        void reset () {
            // Nothing to reset
        }
        void run () {
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
        void run_selected(){
            selected(box,0.1,4);
        }

    private:
        Polytope* box; 
        xyz_f_t selected;
};


#endif // BENCHMARK_H
