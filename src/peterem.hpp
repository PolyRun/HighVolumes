#include <iostream>

extern "C" { // must be included C stlye
#include "volume.h"
}

#ifndef HEADER_PETEREM_HPP
#define HEADER_PETEREM_HPP

void hello(Polytope* p) {
   std::cout << "hello" << std::endl;
   std::cout << "dim: " << p->n << " " << p->m << std::endl;
};

#endif // HEADER_PETEREM_HPP
