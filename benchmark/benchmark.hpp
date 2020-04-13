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

#endif // BENCHMARK_H
