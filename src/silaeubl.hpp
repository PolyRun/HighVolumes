#include <iostream>

extern "C" { // must be included C stlye
#include "random/prng.h"
}

#ifndef HEADER_SILAEUBL_HPP
#define HEADER_SILAEUBL_HPP

void hello() {
    std::cout << "We are goint to print some random numbers:\n" << std::endl;
};

#endif // HEADER_SILAEUBL_HPP
