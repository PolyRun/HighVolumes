#include <iostream>

extern "C" { // must be included C stlye
#include "poly/volume.h"
}

#include "poly/volume_helper.hpp"

#ifndef HEADER_PETEREM_HPP
#define HEADER_PETEREM_HPP

void hello(Polytope* p) {
   std::cout << "hello" << std::endl;
   std::cout << "dim: " << p->n << " " << p->m << std::endl;

   std::cout << *p << std::endl;
   std::cout << p << std::endl;
};

#endif // HEADER_PETEREM_HPP
