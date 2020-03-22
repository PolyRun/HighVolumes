#include <iostream>
#include <cassert>

extern "C" { // must be included C stlye
#include "volume.h"
}


#ifndef HEADER_VOLUME_HELPER_HPP
#define HEADER_VOLUME_HELPER_HPP


std::ostream& operator<<(std::ostream& os, const Polytope& p);
std::ostream& operator<<(std::ostream& os, const Polytope* p);

Polytope* Polytope_new_box(int n, int r);
// allocates / generates cube polytope with radius r


#endif // HEADER_VOLUME_HELPER_HPP



