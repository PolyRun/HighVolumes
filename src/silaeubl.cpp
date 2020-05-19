#include <iostream>
#include "volume/volume_helper.hpp"

extern "C" { // must be included C stlye
#include "random/prng.h"
}

#include "silaeubl.hpp"

int main(int argc, char *argv[]) {
    CLI cli(argc,argv,"benchmark");
    CLIFunctionsVolume cliFun(cli);

    cliFun.preParse();
    if (!cli.parse()) {return -1;}
    cliFun.postParse();

    hello();

    prng_init();
    std::cout << "random doubles:" << std::endl;
    for (int i = 0; i < 32; ++i){
        std::cout <<prng_get_random_double() << std::endl;;
    }
    std::cout << "\nrandom doubles in range -2.0, 6.0:\n";
    for (int i = 0; i < 32; ++i){
        std::cout << " "<< prng_get_random_double_in_range(-2.0, 6.0);
    }
    std::cout << "\nrandom doubles in 0, 1:" << std::endl;
    for (int i = 0; i < 32; ++i){
        std::cout << " " << prng_get_random_double_0_1();
    }
    std::cout << "\nrandom ints:" << std::endl;
    for (int i = 0; i < 32; ++i){
        std::cout << " " << prng_get_random_int();
    }
    std::cout << "\nrandom ints in range -6, 6:" << std::endl;
    for (int i = 0; i < 32; ++i){
        std::cout << " " << prng_get_random_int_in_range(-6, 6);
    }
    std::cout << std::endl;
}

