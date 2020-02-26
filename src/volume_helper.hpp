#include <iostream>

extern "C" { // must be included C stlye
#include "volume.h"
}


#ifndef HEADER_VOLUME_HELPER_HPP
#define HEADER_VOLUME_HELPER_HPP


std::ostream& operator<<(std::ostream& os, const Polytope& p);
std::ostream& operator<<(std::ostream& os, const Polytope* p);




#endif // HEADER_VOLUME_HELPER_HPP



