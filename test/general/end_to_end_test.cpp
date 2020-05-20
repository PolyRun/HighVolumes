#include <iostream>
#include <cassert>

//#include "../../src/volume/volume.h"
#include "../../src/volume/volume_helper.hpp"
//#include "../../src/volume/volume_examples.hpp"
#include "../../src/random/prng.h"


int main(int argc, char **argv) {

    for (int i = 0; i < 10; i++) {
        double num = prng_fast_32_get_random_double_0_1();
        std::cout << num << std::endl;
    }
    
}
