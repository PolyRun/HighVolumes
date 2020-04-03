#include <iostream>
#include <cassert>

extern "C" { // must be included C stlye
#include "volume.h"
}

#include "../util/cli_functions.hpp"

#ifndef HEADER_VOLUME_HELPER_HPP
#define HEADER_VOLUME_HELPER_HPP


std::ostream& operator<<(std::ostream& os, const Polytope& p);
std::ostream& operator<<(std::ostream& os, const Polytope* p);

Polytope* Polytope_new_box(int n, int r);
// allocates / generates cube polytope with radius r


class CLIFunctionsVolume : public CLIFunctions {
public:
   CLIFunctionsVolume(CLI &cli) : CLIFunctions(cli) {
      // please add your functions below.
      // 
      // handle to global variable used to reference current choice of function
      // opt code for cli arg, under which the function is configured (char)
      // string for function name
      // string for default function name
      // map: function name -> function ptr

      add(new CLIF_Option<xyz_f_t>(&xyz_f,'f',"xyz_f","xyz_f1", std::map<std::string, xyz_f_t>{
                                                     {"xyz_f1",xyz_f1},
						     {"xyz_f2",xyz_f2} }));

   }
};



#endif // HEADER_VOLUME_HELPER_HPP



