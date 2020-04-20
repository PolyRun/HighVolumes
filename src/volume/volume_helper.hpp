#include <iostream>
#include <cassert>
#include <vector>
#include <random>
#include <fstream>
#include <glpk.h>
#include <chrono>



extern "C" { // must be included C stlye
#include "volume.h"
}

#include "../util/cli_functions.hpp"
#include "../../polyvest/vol.h"

#ifndef HEADER_VOLUME_HELPER_HPP
#define HEADER_VOLUME_HELPER_HPP


std::ostream& operator<<(std::ostream& os, const Polytope& p);
std::ostream& operator<<(std::ostream& os, const Polytope* p);

Polytope* Polytope_new_box(int n, FT r);
// allocates / generates cube polytope with radius r


class CLIFunctionsVolume : public CLIFunctions {
public:
   CLIFunctionsVolume(CLI &cli) : CLIFunctions(cli) {
      // initialize memory arrays of volume library
      volume_lib_init(200,20);// max_n=200, max_b=20

      // please add your functions below.
      // 
      // handle to global variable used to reference current choice of function
      // opt code for cli arg, under which the function is configured (char)
      // string for function name
      // string for default function name
      // map: function name -> function ptr
      
      claimOpt('f',"Algorithm Functions");
      add(new CLIF_Option<xyz_f_t>(&xyz_f,'f',"xyz_f","xyz_f1", std::map<std::string, xyz_f_t>{
                                                     {"xyz_f1",xyz_f1},
						     {"xyz_f2",xyz_f2} }));

      add(new CLIF_Option<walk_f_t>(&walk_f,'f',"walk_f","walk_ref", std::map<std::string, walk_f_t>{
                                                     {"walk_ref",walk_ref},
						     {"walkCoord_ref",walkCoord_ref} }));

      add(new CLIF_Option<intersectCoord_f_t>(&Polytope_T.intersectCoord,'f',"Polytope_intersectCoord","cached_ref", std::map<std::string, intersectCoord_f_t>{
                                                     {"ref",Polytope_intersectCoord_ref},
						     {"cached_ref",Polytope_intersectCoord_cached_ref} }));

      // number parameters:
      claimOpt('c',"Algorithm Constants");
      add(new CLIF_OptionNumber<int>(&step_size,'c',"step_size","100000", 100, 1e7));
      add(new CLIF_OptionNumber<int>(&walk_size,'c',"walk_size","1", 1, 1e6));
   }
};




/**
 *\brief this function allow creating a random polytope around an ellipsoid given by ell
 * \param ell an n vector standing for a diagonal matrix (note dimension n is implicit in length of ell)
 * \param m the number of constraints
 * \param ret the return polytope
**/
void make_random_poly(const std::vector<double> &ell, int m, Polytope **ret);


/**
 * \brief read in one of the example polytopes defined in polyvest/examples
 * \param filename the relative path from root directory, something like polyvest/examples/simplex_10
 * \param P a *P will be filled with the example polytope from file
 **/
int read_polyvest_p(std::string filename, Polytope **P);


/**
 * \brief convert from Polytope to Polyvest_p
 * \param P the given polytope
 * \param Q is filled with content equivalend to P
 **/
void polyvest_convert(Polytope *P, vol::Polyvest_p *Q);


#endif // HEADER_VOLUME_HELPER_HPP



