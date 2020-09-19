#ifndef HEADER_VOLUME_HELPER_HPP
#define HEADER_VOLUME_HELPER_HPP

#include <iostream>
#include <cassert>
#include <vector>
#include <random>
#include <fstream>
#include <glpk.h>
#include <chrono>

#include <regex>
#include <sstream>


extern "C" { // must be included C stlye
#include "volume.h"
}

#include "../util/cli_functions.hpp"
#include "../util/performance_counter.hpp"
#include "../util/union_find.hpp"
#include "../../polyvest/vol.h"
#include "volume_cost.hpp"
#include "volume_examples.hpp"

std::ostream& operator<<(std::ostream& os, const Polytope& p);
std::ostream& operator<<(std::ostream& os, const Polytope* p);

Polytope* Polytope_new_box(int n, FT r);
PolytopeT* PolytopeT_new_box(int n, FT r);
// allocates / generates cube polytope with radius r


/**
 * sets up all the options for function choice.
 * Run this right after you have created the clif
 **/
void initVolumeFunctions(CLIFunctions &clif);


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
int read_vinci(string filename, Polytope **P, FT *vol);


/**
 * \brief convert from Polytope to Polyvest_p
 * \param P the given polytope
 * \param Q is filled with content equivalend to P
 **/
void polyvest_convert(Polytope *P, vol::Polyvest_p *Q);


// tries to shuffle rows such, that there are the most pairs of non-zeros
// pair: two non-zero entries in the same column, that lie one row apart
// this should lead to better access and vectorization
Polytope* optimize_polytope(Polytope *p);


#endif // HEADER_VOLUME_HELPER_HPP



