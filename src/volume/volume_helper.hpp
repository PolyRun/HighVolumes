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


class CLIFunctionsVolume : public CLIFunctions {
public:
   CLIFunctionsVolume(CLI &cli) : CLIFunctions(cli) {
      // initialize memory arrays of volume library
       // note, cross_13 has 8000 constraints so let's set max_m to that
       volume_lib_init(200,20000, 20);// max_n=200, max_m=20000, max_b=20
      
      // add performance counter cost functions
      pc_stack().add((void*)dotProduct_ref, new PC_Cost_Wrapper<dotProduct_cost_f>(dotProduct_cost_ref,"dotProduct_ref"));
      pc_stack().add((void*)dotProduct_2acc, new PC_Cost_Wrapper<dotProduct_cost_f>(dotProduct_cost_ref,"dotProduct_2acc"));
      pc_stack().add((void*)dotProduct_auto1, new PC_Cost_Wrapper<dotProduct_cost_f>(dotProduct_cost_ref,"dotProduct_auto1"));
      pc_stack().add((void*)dotProduct_auto2, new PC_Cost_Wrapper<dotProduct_cost_f>(dotProduct_cost_ref,"dotProduct_auto2"));
      pc_stack().add((void*)dotProduct_vec1, new PC_Cost_Wrapper<dotProduct_cost_f>(dotProduct_cost_ref,"dotProduct_vec1"));
      
      pc_stack().add((void*)squaredNorm, new PC_Cost_Wrapper<squaredNorm_cost_f>(squaredNorm_cost_ref,"squaredNorm_ref"));
      
      pc_stack().add((void*)Ball_intersect_ref, new PC_Cost_Wrapper<Ball_intersect_cost_f>(Ball_intersect_cost_ref,"Ball_intersect_ref"));
      pc_stack().add((void*)Ball_intersectCoord_ref, new PC_Cost_Wrapper<Ball_intersectCoord_cost_f>(Ball_intersectCoord_cost_ref,"Ball_intersectCoord_ref"));
      
      pc_stack().add((void*)Ball_intersectCoord_cached_ref, new PC_Cost_Wrapper<Ball_intersectCoord_cached_cost_f>(Ball_intersectCoord_cached_cost_ref,"Ball_intersectCoord_cached_ref"));
      pc_stack().add((void*)Ball_intersectCoord_cached4, new PC_Cost_Wrapper<Ball_intersectCoord_cached_cost_f>(Ball_intersectCoord_cached4_cost_ref,"Ball_intersectCoord_cached4"));
      pc_stack().add((void*)Ball_intersectCoord_cached8, new PC_Cost_Wrapper<Ball_intersectCoord_cached_cost_f>(Ball_intersectCoord_cached8_cost_ref,"Ball_intersectCoord_cached8"));
      
      // Polytope
      pc_stack().add((void*)Polytope_intersect_ref, new PC_Cost_Wrapper<intersect_cost_f>(Polytope_intersect_cost_ref,"Polytope_intersect_ref"));
      pc_stack().add((void*)Polytope_intersectCoord_ref, new PC_Cost_Wrapper<intersectCoord_cost_f>(Polytope_intersectCoord_cost_ref,"Polytope_intersectCoord_ref"));
      pc_stack().add((void*)Polytope_intersectCoord_cached_ref, new PC_Cost_Wrapper<intersectCoord_cost_f>(Polytope_intersectCoord_cached_cost_ref,"Polytope_intersectCoord_cached_ref"));
      pc_stack().add((void*)Polytope_cacheUpdateCoord_ref, new PC_Cost_Wrapper<cacheUpdateCoord_cost_f>(Polytope_cacheUpdateCoord_cost_ref,"Polytope_cacheUpdateCoord_ref")); 
      pc_stack().add((void*)Polytope_cacheReset_ref, new PC_Cost_Wrapper<cacheReset_cost_f>(Polytope_cacheReset_cost_ref,"Polytope_cacheReset_ref"));
      
      // PolytopeT
      pc_stack().add((void*)PolytopeT_intersect_ref, new PC_Cost_Wrapper<intersect_cost_f>(PolytopeT_intersect_cost_ref,"PolytopeT_intersect_ref"));
      
      pc_stack().add((void*)PolytopeT_intersectCoord_cached_b_vec, new PC_Cost_Wrapper<intersectCoord_cost_f>(PolytopeT_intersectCoord_cached_b_cost_ref,"PolytopeT_intersectCoord_cached_b_vec"));
      pc_stack().add((void*)PolytopeT_intersectCoord_cached_b_vec_inl, new PC_Cost_Wrapper<intersectCoord_cost_f>(PolytopeT_intersectCoord_cached_b_cost_ref,"PolytopeT_intersectCoord_cached_b_vec_inl"));
      pc_stack().add((void*)PolytopeT_intersectCoord_cached_b_vec2, new PC_Cost_Wrapper<intersectCoord_cost_f>(PolytopeT_intersectCoord_cached_b_cost_ref,"PolytopeT_intersectCoord_cached_b_vec2"));
      pc_stack().add((void*)PolytopeT_intersectCoord_ref, new PC_Cost_Wrapper<intersectCoord_cost_f>(PolytopeT_intersectCoord_cost_ref,"PolytopeT_intersectCoord_ref"));
      
      pc_stack().add((void*)PolytopeT_intersectCoord_cached_ref, new PC_Cost_Wrapper<intersectCoord_cost_f>(PolytopeT_intersectCoord_cached_cost_ref,"PolytopeT_intersectCoord_cached_ref"));
      pc_stack().add((void*)PolytopeT_intersectCoord_cached_nc1, new PC_Cost_Wrapper<intersectCoord_cost_f>(PolytopeT_intersectCoord_cached_cost_ref,"PolytopeT_intersectCoord_cached_nc1"));
      pc_stack().add((void*)PolytopeT_intersectCoord_cached_b_ref, new PC_Cost_Wrapper<intersectCoord_cost_f>(PolytopeT_intersectCoord_cached_b_cost_ref,"PolytopeT_intersectCoord_cached_b_ref"));
      pc_stack().add((void*)PolytopeT_intersectCoord_vectorized, new PC_Cost_Wrapper<intersectCoord_cost_f>(PolytopeT_intersectCoord_cached_cost_ref,"PolytopeT_intersectCoord_vectorized"));
      
      pc_stack().add((void*)PolytopeT_intersectCoord_cached_b_inv_ref, new PC_Cost_Wrapper<intersectCoord_cost_f>(PolytopeT_intersectCoord_cached_b_cost_ref,"PolytopeT_intersectCoord_cached_b_inv_ref"));
      pc_stack().add((void*)PolytopeT_intersectCoord_cached_b_inv_vec, new PC_Cost_Wrapper<intersectCoord_cost_f>(PolytopeT_intersectCoord_cached_b_cost_ref,"PolytopeT_intersectCoord_cached_b_inv_vec"));
      pc_stack().add((void*)PolytopeT_intersectCoord_cached_b_inv_vec, new PC_Cost_Wrapper<intersectCoord_cost_f>(PolytopeT_intersectCoord_cached_b_cost_ref,"PolytopeT_intersectCoord_cached_b_inv_vec_inl"));
      
      pc_stack().add((void*)PolytopeT_cacheUpdateCoord_ref, new PC_Cost_Wrapper<cacheUpdateCoord_cost_f>(PolytopeT_cacheUpdateCoord_cost_ref,"PolytopeT_cacheUpdateCoord_ref"));
      pc_stack().add((void*)PolytopeT_cacheUpdateCoord_b_ref, new PC_Cost_Wrapper<cacheUpdateCoord_cost_f>(PolytopeT_cacheUpdateCoord_b_cost_ref,"PolytopeT_cacheUpdateCoord_b_ref"));
      pc_stack().add((void*)PolytopeT_cacheUpdateCoord_b_vec, new PC_Cost_Wrapper<cacheUpdateCoord_cost_f>(PolytopeT_cacheUpdateCoord_b_cost_vec,"PolytopeT_cacheUpdateCoord_b_vec"));
      
      pc_stack().add((void*)PolytopeT_cacheReset_ref, new PC_Cost_Wrapper<cacheReset_cost_f>(PolytopeT_cacheReset_cost_ref,"PolytopeT_cacheReset_ref"));
      pc_stack().add((void*)PolytopeT_cacheReset_b_ref, new PC_Cost_Wrapper<cacheReset_cost_f>(PolytopeT_cacheReset_b_cost_ref,"PolytopeT_cacheReset_b_ref"));
      pc_stack().add((void*)PolytopeT_cacheReset_b_vec, new PC_Cost_Wrapper<cacheReset_cost_f>(PolytopeT_cacheReset_b_cost_vec,"PolytopeT_cacheReset_b_vec"));
      
      pc_stack().add((void*)PolytopeT_intersectCoord4_ref, new PC_Cost_Wrapper<intersectCoord_cost_f>(PolytopeT_intersectCoord4_cost_ref,"PolytopeT_intersectCoord4_ref"));
      pc_stack().add((void*)PolytopeT_intersectCoord8_ref, new PC_Cost_Wrapper<intersectCoord_cost_f>(PolytopeT_intersectCoord8_cost_ref,"PolytopeT_intersectCoord8_ref"));
      pc_stack().add((void*)PolytopeT_cacheUpdateCoord4_ref, new PC_Cost_Wrapper<cacheUpdateCoord_cost_f>(PolytopeT_cacheUpdateCoord4_cost_ref,"PolytopeT_cacheUpdateCoord4_ref"));
      pc_stack().add((void*)PolytopeT_cacheUpdateCoord8_ref, new PC_Cost_Wrapper<cacheUpdateCoord_cost_f>(PolytopeT_cacheUpdateCoord8_cost_ref,"PolytopeT_cacheUpdateCoord8_ref"));
      pc_stack().add((void*)PolytopeT_cacheReset4_ref, new PC_Cost_Wrapper<cacheReset_cost_f>(PolytopeT_cacheReset4_cost_ref,"PolytopeT_cacheReset4_ref"));
      pc_stack().add((void*)PolytopeT_cacheReset8_ref, new PC_Cost_Wrapper<cacheReset_cost_f>(PolytopeT_cacheReset8_cost_ref,"PolytopeT_cacheReset8_ref"));
 
      // Ellipsoid
      pc_stack().add((void*)Ellipsoid_intersect_ref, new PC_Cost_Wrapper<intersect_cost_f>(Ellipsoid_intersect_cost_ref,"Ellipsoid_intersect_ref"));
      
      pc_stack().add((void*)Ellipsoid_intersectCoord_ref, new PC_Cost_Wrapper<intersectCoord_cost_f>(Ellipsoid_intersectCoord_cost_ref,"Ellipsoid_intersectCoord_ref"));
      pc_stack().add((void*)Ellipsoid_intersectCoord_cached_ref, new PC_Cost_Wrapper<intersectCoord_cost_f>(Ellipsoid_intersectCoord_cached_cost_ref,"Ellipsoid_intersectCoord_cached_ref"));
      pc_stack().add((void*)Ellipsoid_intersectCoord_cached_reord, new PC_Cost_Wrapper<intersectCoord_cost_f>(Ellipsoid_intersectCoord_cached_cost_ref,"Ellipsoid_intersectCoord_cached_reord"));
      
      pc_stack().add((void*)Ellipsoid_intersectCoord_cached_reord2, new PC_Cost_Wrapper<intersectCoord_cost_f>(Ellipsoid_intersectCoord_cached_cost_ref,"Ellipsoid_intersectCoord_cached_reord2"));
      pc_stack().add((void*)Ellipsoid_intersectCoord_cached_reord3, new PC_Cost_Wrapper<intersectCoord_cost_f>(Ellipsoid_intersectCoord_cached_cost_ref,"Ellipsoid_intersectCoord_cached_reord3"));
      pc_stack().add((void*)Ellipsoid_intersectCoord_cached_reord_fma, new PC_Cost_Wrapper<intersectCoord_cost_f>(Ellipsoid_intersectCoord_cached_cost_reord_fma,"Ellipsoid_intersectCoord_cached_reord_fma"));
      pc_stack().add((void*)Ellipsoid_cacheUpdateCoord_ref, new PC_Cost_Wrapper<cacheUpdateCoord_cost_f>(Ellipsoid_cacheUpdateCoord_cost_ref,"Ellipsoid_cacheUpdateCoord_ref"));
      pc_stack().add((void*)Ellipsoid_cacheUpdateCoord_c, new PC_Cost_Wrapper<cacheUpdateCoord_cost_f>(Ellipsoid_cacheUpdateCoord_cost_ref,"Ellipsoid_cacheUpdateCoord_c"));
      pc_stack().add((void*)Ellipsoid_cacheUpdateCoord_fma, new PC_Cost_Wrapper<cacheUpdateCoord_cost_f>(Ellipsoid_cacheUpdateCoord_cost_ref,"Ellipsoid_cacheUpdateCoord_fma"));
      
      pc_stack().add((void*)Ellipsoid_cacheUpdateCoord_vec, new PC_Cost_Wrapper<cacheUpdateCoord_cost_f>(Ellipsoid_cacheUpdateCoord_cost_ref,"Ellipsoid_cacheUpdateCoord_vec"));
      pc_stack().add((void*)Ellipsoid_cacheUpdateCoord_vec_u2, new PC_Cost_Wrapper<cacheUpdateCoord_cost_f>(Ellipsoid_cacheUpdateCoord_cost_ref,"Ellipsoid_cacheUpdateCoord_vec_u2"));
      pc_stack().add((void*)Ellipsoid_cacheUpdateCoord_vec_u4, new PC_Cost_Wrapper<cacheUpdateCoord_cost_f>(Ellipsoid_cacheUpdateCoord_cost_ref,"Ellipsoid_cacheUpdateCoord_vec_u4"));
      
      pc_stack().add((void*)Ellipsoid_cacheUpdateCoord_vec2, new PC_Cost_Wrapper<cacheUpdateCoord_cost_f>(Ellipsoid_cacheUpdateCoord_cost_ref,"Ellipsoid_cacheUpdateCoord_vec2"));
      pc_stack().add((void*)Ellipsoid_cacheUpdateCoord_vec2_u2, new PC_Cost_Wrapper<cacheUpdateCoord_cost_f>(Ellipsoid_cacheUpdateCoord_cost_ref,"Ellipsoid_cacheUpdateCoord_vec2_u2"));
      pc_stack().add((void*)Ellipsoid_cacheUpdateCoord_vec2_u4, new PC_Cost_Wrapper<cacheUpdateCoord_cost_f>(Ellipsoid_cacheUpdateCoord_cost_ref,"Ellipsoid_cacheUpdateCoord_vec2_u4"));
      
      pc_stack().add((void*)Ellipsoid_cacheReset_ref, new PC_Cost_Wrapper<cacheReset_cost_f>(Ellipsoid_cacheReset_cost_ref,"Ellipsoid_cacheReset_ref"));
      pc_stack().add((void*)Ellipsoid_cacheReset_reord, new PC_Cost_Wrapper<cacheReset_cost_f>(Ellipsoid_cacheReset_cost_ref,"Ellipsoid_cacheReset_reord"));
      pc_stack().add((void*)Ellipsoid_cacheReset_fma, new PC_Cost_Wrapper<cacheReset_cost_f>(Ellipsoid_cacheReset_cost_ref,"Ellipsoid_cacheReset_fma"));

      // PolytopeCSC
      pc_stack().add((void *) PolytopeCSC_intersect_ref, new PC_Cost_Wrapper<intersect_cost_f>(PolytopeCSC_intersect_cost_ref, "PolytopeCSC_intersect_ref"));
      pc_stack().add((void *) PolytopeCSC_mvm, new PC_Cost_Wrapper<mvm_cost_f>(PolytopeCSC_mvm_cost, "PolytopeCSC_mvm"));
      
      pc_stack().add((void *) PolytopeCSC_intersectCoord_ref, new PC_Cost_Wrapper<intersectCoord_cost_f>(PolytopeCSC_intersectCoord_cost_ref, "PolytopeCSC_intersectCoord_ref"));
      pc_stack().add((void *) PolytopeCSC_intersectCoord_cached_ref, new PC_Cost_Wrapper<intersectCoord_cost_f>(PolytopeCSC_intersectCoord_cached_cost_ref, "PolytopeCSC_intersectCoord_cached_ref"));
      pc_stack().add((void *) PolytopeCSC_intersectCoord_cached_withb, new PC_Cost_Wrapper<intersectCoord_cost_f>(PolytopeCSC_intersectCoord_cached_cost_withb, "PolytopeCSC_intersectCoord_cached_withb"));
      pc_stack().add((void *) PolytopeCSC_intersectCoord_cached_vec, new PC_Cost_Wrapper<intersectCoord_cost_f>(PolytopeCSC_intersectCoord_cached_cost_vec, "PolytopeCSC_intersectCoord_cached_vec"));
      pc_stack().add((void *) PolytopeCSC_intersectCoord_cached_vec_nogather, new PC_Cost_Wrapper<intersectCoord_cost_f>(PolytopeCSC_intersectCoord_cached_cost_vec, "PolytopeCSC_intersectCoord_cached_vec_nogather"));
      pc_stack().add((void *) PolytopeCSC_intersectCoord_cached_vec_nan, new PC_Cost_Wrapper<intersectCoord_cost_f>(PolytopeCSC_intersectCoord_cached_cost_vec, "PolytopeCSC_intersectCoord_cached_vec_nan"));
      pc_stack().add((void *) PolytopeCSC_intersectCoord_cached_vec_nan_inv, new PC_Cost_Wrapper<intersectCoord_cost_f>(PolytopeCSC_intersectCoord_cached_cost_vec, "PolytopeCSC_intersectCoord_cached_vec_nan_inv"));
      pc_stack().add((void *) PolytopeCSC_intersectCoord_cached_vec_inline, new PC_Cost_Wrapper<intersectCoord_cost_f>(PolytopeCSC_intersectCoord_cached_cost_vec, "PolytopeCSC_intersectCoord_cached_vec_inline"));
      pc_stack().add((void *) PolytopeCSC_intersectCoord_cached_vec_inline_2accs, new PC_Cost_Wrapper<intersectCoord_cost_f>(PolytopeCSC_intersectCoord_cached_cost_vec, "PolytopeCSC_intersectCoord_cached_cost_vec_inline_2accs"));
      
      pc_stack().add((void *) PolytopeCSC_cacheReset_ref, new PC_Cost_Wrapper<cacheReset_cost_f>(PolytopeCSC_cacheReset_cost_ref, "PolytopeCSC_cacheReset_ref"));
      pc_stack().add((void *) PolytopeCSC_cacheReset_withb, new PC_Cost_Wrapper<cacheReset_cost_f>(PolytopeCSC_cacheReset_cost_withb, "PolytopeCSC_cacheReset_withb"));
      pc_stack().add((void *) PolytopeCSC_cacheReset_fma, new PC_Cost_Wrapper<cacheReset_cost_f>(PolytopeCSC_cacheReset_cost_withb, "PolytopeCSC_cacheReset_fma"));
      pc_stack().add((void *) PolytopeCSC_cacheReset_vec, new PC_Cost_Wrapper<cacheReset_cost_f>(PolytopeCSC_cacheReset_cost_withb, "PolytopeCSC_cacheReset_vec"));
      
      pc_stack().add((void *) PolytopeCSC_cacheUpdateCoord_ref, new PC_Cost_Wrapper<cacheUpdateCoord_cost_f>(PolytopeCSC_cacheUpdateCoord_cost_ref, "PolytopeCSC_cacheUpdateCoord_ref"));
      pc_stack().add((void *) PolytopeCSC_cacheUpdateCoord_withb, new PC_Cost_Wrapper<cacheUpdateCoord_cost_f>(PolytopeCSC_cacheUpdateCoord_cost_withb, "PolytopeCSC_cacheUpdateCoord_withb"));
      pc_stack().add((void *) PolytopeCSC_cacheUpdateCoord_fma, new PC_Cost_Wrapper<cacheUpdateCoord_cost_f>(PolytopeCSC_cacheUpdateCoord_cost_withb, "PolytopeCSC_cacheUpdateCoord_fma"));
      pc_stack().add((void *) PolytopeCSC_cacheUpdateCoord_vec, new PC_Cost_Wrapper<cacheUpdateCoord_cost_f>(PolytopeCSC_cacheUpdateCoord_cost_withb, "PolytopeCSC_cacheUpdateCoord_vec"));

      pc_stack().add((void*)PolytopeCSC_intersectCoord4_ref, new PC_Cost_Wrapper<intersectCoord_cost_f>(PolytopeCSC_intersectCoord4_cost_ref,"PolytopeCSC_intersectCoord4_ref"));
      pc_stack().add((void*)PolytopeCSC_intersectCoord8_ref, new PC_Cost_Wrapper<intersectCoord_cost_f>(PolytopeCSC_intersectCoord8_cost_ref,"PolytopeCSC_intersectCoord8_ref"));
      pc_stack().add((void*)PolytopeCSC_cacheUpdateCoord4_ref, new PC_Cost_Wrapper<cacheUpdateCoord_cost_f>(PolytopeCSC_cacheUpdateCoord4_cost_ref,"PolytopeCSC_cacheUpdateCoord4_ref"));
      pc_stack().add((void*)PolytopeCSC_cacheUpdateCoord8_ref, new PC_Cost_Wrapper<cacheUpdateCoord_cost_f>(PolytopeCSC_cacheUpdateCoord8_cost_ref,"PolytopeCSC_cacheUpdateCoord8_ref"));
      pc_stack().add((void*)PolytopeCSC_cacheReset4_ref, new PC_Cost_Wrapper<cacheReset_cost_f>(PolytopeCSC_cacheReset4_cost_ref,"PolytopeCSC_cacheReset4_ref"));
      pc_stack().add((void*)PolytopeCSC_cacheReset8_ref, new PC_Cost_Wrapper<cacheReset_cost_f>(PolytopeCSC_cacheReset8_cost_ref,"PolytopeCSC_cacheReset8_ref"));

      pc_stack().add((void *) PolytopeCSC_intersectCoord_cached_vec_onlyread, new PC_Cost_Wrapper<intersectCoord_cost_f>(PolytopeCSC_intersectCoord_cached_cost_vec, "PolytopeCSC_intersectCoord_cached_vec_onlyread"));
       
      // PolytopeJIT
      pc_stack().add((void*)PolytopeJIT_intersect_ref, new PC_Cost_Wrapper<intersect_cost_f>(PolytopeJIT_intersect_cost_ref,"PolytopeJIT_intersect_ref"));
      pc_stack().add((void*)PolytopeJIT_intersectCoord_ref, new PC_Cost_Wrapper<intersectCoord_cost_f>(PolytopeJIT_intersectCoord_cost_ref,"PolytopeJIT_intersectCoord_ref"));
      pc_stack().add((void*)PolytopeJIT_cacheUpdateCoord_ref, new PC_Cost_Wrapper<cacheUpdateCoord_cost_f>(PolytopeJIT_cacheUpdateCoord_cost_ref,"PolytopeJIT_cacheUpdateCoord_ref"));
      pc_stack().add((void*)PolytopeJIT_cacheReset_ref, new PC_Cost_Wrapper<cacheReset_cost_f>(PolytopeJIT_cacheReset_cost_ref,"PolytopeJIT_cacheReset_ref"));

      pc_stack().add((void*)PolytopeJIT_intersectCoord4_ref, new PC_Cost_Wrapper<intersectCoord_cost_f>(PolytopeJIT_intersectCoord4_cost_ref,"PolytopeJIT_intersectCoord4_ref"));
      pc_stack().add((void*)PolytopeJIT_intersectCoord8_ref, new PC_Cost_Wrapper<intersectCoord_cost_f>(PolytopeJIT_intersectCoord8_cost_ref,"PolytopeJIT_intersectCoord8_ref"));
      pc_stack().add((void*)PolytopeJIT_cacheUpdateCoord4_ref, new PC_Cost_Wrapper<cacheUpdateCoord_cost_f>(PolytopeJIT_cacheUpdateCoord4_cost_ref,"PolytopeJIT_cacheUpdateCoord4_ref"));
      pc_stack().add((void*)PolytopeJIT_cacheUpdateCoord8_ref, new PC_Cost_Wrapper<cacheUpdateCoord_cost_f>(PolytopeJIT_cacheUpdateCoord8_cost_ref,"PolytopeJIT_cacheUpdateCoord8_ref"));
      pc_stack().add((void*)PolytopeJIT_cacheReset4_ref, new PC_Cost_Wrapper<cacheReset_cost_f>(PolytopeJIT_cacheReset4_cost_ref,"PolytopeJIT_cacheReset4_ref"));
      pc_stack().add((void*)PolytopeJIT_cacheReset8_ref, new PC_Cost_Wrapper<cacheReset_cost_f>(PolytopeJIT_cacheReset8_cost_ref,"PolytopeJIT_cacheReset8_ref"));

      // volume
      pc_stack().add((void*)volume_ref, new PC_Cost_Wrapper<volume_cost_f>(volume_cost_ref,"volume_ref"));
      pc_stack().add((void*)volume_coord_single, new PC_Cost_Wrapper<volume_cost_f>(volume_coord_1_cost_ref,"volume_coord_1"));
      pc_stack().add((void*)volume_coord_4,      new PC_Cost_Wrapper<volume_cost_f>(volume_coord_4_cost_ref,"volume_coord_4"));
      pc_stack().add((void*)volume_coord_8,      new PC_Cost_Wrapper<volume_cost_f>(volume_coord_8_cost_ref,"volume_coord_8"));
      
      pc_stack().add((void*)walk_ref, new PC_Cost_Wrapper<walk_cost_f>(walk_cost_ref,"walk_ref"));
      pc_stack().add((void*)walkCoord_ref, new PC_Cost_Wrapper<walk_cost_f>(walkCoord_cost_ref,"walkCoord_ref"));
      pc_stack().add((void*)walkCoord_coord_single, new PC_Cost_Wrapper<walk_cost_f>(walkCoord_coord_1_cost_ref,"walkCoord_coord_1"));
      pc_stack().add((void*)walkCoord_coord_4,      new PC_Cost_Wrapper<walk_cost_f>(walkCoord_coord_4_cost_ref,"walkCoord_coord_4"));
      pc_stack().add((void*)walkCoord_coord_8,      new PC_Cost_Wrapper<walk_cost_f>(walkCoord_coord_8_cost_ref,"walkCoord_coord_8"));

      // Randomness
      pc_stack().add((void *)prng_get_random_int, new PC_Cost_Wrapper<random_int_cost_f>(Random_int_cost_ref, "Random int"));
      pc_stack().add((void *)prng_get_random_int_in_range, new PC_Cost_Wrapper<random_int_in_range_cost_f>(Random_int_in_range_cost_ref, "Random int_in_range"));
      pc_stack().add((void *)prng_get_random_double_in_range, new PC_Cost_Wrapper<random_double_in_range_cost_f>(Random_double_in_range_cost_ref, "Random double_in_range"));
      pc_stack().add((void *)prng_get_random_double_0_1, new PC_Cost_Wrapper<random_double_0_1_cost_f>(Random_double_0_1_cost_ref, "Random double_0_1"));
      pc_stack().add((void *)prng_get_random_double_normal, new PC_Cost_Wrapper<random_double_normal_cost_f>(Random_double_normal_cost_ref, "Random double_normal"));
      
      pc_stack().add((void *)std_rand256d, new PC_Cost_Wrapper<rand256d_cost_f_t>(sr_rand256d_cost_ref, "Random 4 doubles (std)"));
      pc_stack().add((void *)sr_rand256d, new PC_Cost_Wrapper<rand256d_cost_f_t>(sr_rand256d_cost_ref, "Random 4 doubles (sr)"));
      
      

      // please add your functions below.
      // 
      // handle to global variable used to reference current choice of function
      // opt code for cli arg, under which the function is configured (char)
      // string for function name
      // string for default function name
      // map: function name -> function ptr
      
      claimOpt('f',"Algorithm Functions");
      add(new CLIF_Option<xyz_f_t>(&xyz_f,'f',"xyz_f","xyz_f1", {
                                                     {"xyz_f1", {xyz_f1, "test 1"}},
						     {"xyz_f2", {xyz_f2, "test 2"}} }));
      add(new CLIF_Option<dotProduct_f_t>(&dotProduct,'f',"dotProduct","2acc", {
                                                     {"ref",  {dotProduct_ref, "ref"}},
                                                     {"2acc", {dotProduct_2acc, "2 accumulators"}},
                                                     {"auto1",{dotProduct_auto1,"auto gen 1"}},
                                                     {"auto2",{dotProduct_auto2,"auto gen 2"}},
						     {"vec1", {dotProduct_vec1,"vectorize"}} }));
      
      add(new CLIF_Option<squaredNorm_f_t>(&squaredNorm,'f',"squaredNorm","ref", {
						     {"ref", {squaredNorm_ref, "ref"}} }));

      add(new CLIF_TrippleOption<squaredNorm_cached_f_t,squaredNorm_cached_reset_f_t,squaredNorm_cached_update_f_t>(
		 &squaredNorm_cached,&squaredNorm_cached_reset,&squaredNorm_cached_update,
		 'f',"squaredNorm_cached","ref", {
                        {"no",        {{squaredNorm_cached_ref,  {squaredNorm_cached_reset_ref, squaredNorm_cached_update_ref }}, "no cacheing (pre-ref)"}},
                        {"ref",       {{squaredNorm_cached_refc, {squaredNorm_cached_reset_refc,squaredNorm_cached_update_refc}}, "cacheing (ref)"}},
		       	}));

      add(new CLIF_DoubleOption<shell_cache_init_f_t, shell_idx_f_t>
          (
           &shell_cache_init, &shell_idx,
           'f', "log_f", "nocache",
           {
            {"nocache", {{shell_cache_init_nocache, shell_idx_nocache}, "no caching (ordinary logs)"}},
            {"ref", {{shell_cache_init_ref, shell_idx_ref}, "linear search in cache"}},
            {"binary", {{shell_cache_init_ref, shell_idx_binary}, "binary search in cache"}},
           }));

                                                                     


      add(new CLIF_DoubleOption<walk_f_t,volume_f_t>(&walk_f,&volume,'f',"walk_f","walkCoord_ref", {
                                                     {"walk_ref",          {{walk_ref,               volume_ref},          "random direction walk (ref)"}},
						     {"walkCoord_ref",     {{walkCoord_ref,          volume_ref},          "coordinate walk (ref)"}},
						     {"walkCoord_1",       {{walkCoord_coord_single, volume_coord_single}, "coordinate walk, cached squaredNorm"}},
						     {"walkCoord_4",       {{walkCoord_coord_4,      volume_coord_4},      "coordinate walk, cached squaredNorm, 4-way parallel x"}},
						     {"walkCoord_8",       {{walkCoord_coord_8,      volume_coord_8},      "coordinate walk, cached squaredNorm, 8-way parallel x"}},
						     }));


      add(new CLIF_Option<intersectCoord_f_t>(&Polytope_T.intersectCoord,'f',"Polytope_intersectCoord","cached_ref", {
                                                     {"ref",        {Polytope_intersectCoord_ref, "no cache (ref)"}},
						     {"cached_ref", {Polytope_intersectCoord_cached_ref, "with cache (ref)"}} }));

      add(new CLIF_TrippleOption<intersectCoord_f_t,cacheReset_f_t,cacheUpdateCoord_f_t>(
		 &PolytopeT_T.intersectCoord,&PolytopeT_T.cacheReset,&PolytopeT_T.cacheUpdateCoord,
		 'f',"PolytopeT_intersectCoord","cached_b_ref", {
                        {"ref",        {{PolytopeT_intersectCoord_ref, {PolytopeT_cacheReset_ref,PolytopeT_cacheUpdateCoord_ref}}, "no cache (ref)"}},
                        {"cached_ref",        {{PolytopeT_intersectCoord_cached_ref, {PolytopeT_cacheReset_ref,PolytopeT_cacheUpdateCoord_ref}}, "with cache (ref)"}},
                        {"cached_nc1",        {{PolytopeT_intersectCoord_cached_nc1, {PolytopeT_cacheReset_ref,PolytopeT_cacheUpdateCoord_ref}}, "with cache, no condition - failed though"}},
                        {"cached_b_vec",        {{PolytopeT_intersectCoord_cached_b_vec, {PolytopeT_cacheReset_b_vec,PolytopeT_cacheUpdateCoord_b_vec}}, "with cache, b and vectorized"}},
                        {"cached_b_vec2",        {{PolytopeT_intersectCoord_cached_b_vec2, {PolytopeT_cacheReset_b_vec,PolytopeT_cacheUpdateCoord_b_vec}}, "with cache, b and vectorized version 2"}},
                        {"cached_b_vec_inl",        {{PolytopeT_intersectCoord_cached_b_vec_inl, {PolytopeT_cacheReset_b_vec,PolytopeT_cacheUpdateCoord_b_vec}}, "with cache, b and vectorized loop unrolling"}},
                        {"cached_b_ref",        {{PolytopeT_intersectCoord_cached_b_ref, {PolytopeT_cacheReset_b_ref,PolytopeT_cacheUpdateCoord_b_ref}}, "with cache, b in cache (ref)"}},
                        {"cached_vectorized", {{PolytopeT_intersectCoord_vectorized, {PolytopeT_cacheReset_ref,PolytopeT_cacheUpdateCoord_ref}}, "with cache and vectorized"}},
                        {"cached_b_inv_ref", {{PolytopeT_intersectCoord_cached_b_inv_ref, {PolytopeT_cacheReset_b_ref,PolytopeT_cacheUpdateCoord_b_ref}}, "with cache b, and inv used for intersect (ref)"}},
                        {"cached_b_inv_vec", {{PolytopeT_intersectCoord_cached_b_inv_vec, {PolytopeT_cacheReset_b_vec,PolytopeT_cacheUpdateCoord_b_vec}}, "with cache b, and inv used for intersect, vectorized"}},
                        {"cached_b_inv_vec_inl", {{PolytopeT_intersectCoord_cached_b_inv_vec_inl, {PolytopeT_cacheReset_b_vec,PolytopeT_cacheUpdateCoord_b_vec}}, "with cache b, and inv used for intersect, vectorized and 2x loop unrolling"}},
		       	}));

      add(new CLIF_TrippleOption<intersectCoord_f_t, cacheReset_f_t, cacheUpdateCoord_f_t>
          (&PolytopeCSC_T.intersectCoord, &PolytopeCSC_T.cacheReset, &PolytopeCSC_T.cacheUpdateCoord,
           'f', "PolytopeCSC_intersectCoord", "cached_b_ref",
           {
            {"ref", {{PolytopeCSC_intersectCoord_ref, {PolytopeCSC_cacheReset_ref, PolytopeCSC_cacheUpdateCoord_ref}}, "no cche (ref)"}},
            {"cached_ref", {{PolytopeCSC_intersectCoord_cached_ref, {PolytopeCSC_cacheReset_ref, PolytopeCSC_cacheUpdateCoord_ref}}, "with cache (ref)"}},
            {"cached_b_ref", {{PolytopeCSC_intersectCoord_cached_withb, {PolytopeCSC_cacheReset_withb, PolytopeCSC_cacheUpdateCoord_withb}}, "with cache, b in cache"}},
            {"cached_b_vec", {{PolytopeCSC_intersectCoord_cached_vec, {PolytopeCSC_cacheReset_fma, PolytopeCSC_cacheUpdateCoord_fma}}, "vectorized and fma, with cache, b in cache"}},
            {"cached_b_vec_nogather", {{PolytopeCSC_intersectCoord_cached_vec_nogather, {PolytopeCSC_cacheReset_fma, PolytopeCSC_cacheUpdateCoord_fma}}, "vectorized (without gather instr) and fma, with cache, b in cache"}},
            {"cached_b_vec_nan", {{PolytopeCSC_intersectCoord_cached_vec_nan, {PolytopeCSC_cacheReset_fma, PolytopeCSC_cacheUpdateCoord_fma}}, "vectorized (without gather instr, use nan) and fma, with cache, b in cache"}},
            {"cached_b_vec_vec_nan", {{PolytopeCSC_intersectCoord_cached_vec_nan, {PolytopeCSC_cacheReset_vec, PolytopeCSC_cacheUpdateCoord_vec}}, "vectorized (without gather instr, use nan) and vector fma in cache update and reset, with cache, b in cache"}},
            {"cached_b_vec_nan_inv", {{PolytopeCSC_intersectCoord_cached_vec_nan_inv, {PolytopeCSC_cacheReset_fma, PolytopeCSC_cacheUpdateCoord_fma}}, "vectorized (without gather instr, use nan) and fma and store Ainv, with cache, b in cache"}},
            {"cached_b_vec_inl", {{PolytopeCSC_intersectCoord_cached_vec_inline, {PolytopeCSC_cacheReset_fma, PolytopeCSC_cacheUpdateCoord_fma}}, "vectorized inlined, with cache, b in cache"}},
            {"cached_b_vec_inl_2accs", {{PolytopeCSC_intersectCoord_cached_vec_inline_2accs, {PolytopeCSC_cacheReset_fma, PolytopeCSC_cacheUpdateCoord_fma}}, "vectorized inlined 2 accs, with cache, b in cache"}},
            // onlyread is meant for testing io bound, doesn't compute the right thing
            {"cached_vec_onlyread", {{PolytopeCSC_intersectCoord_cached_vec_onlyread, {PolytopeCSC_cacheReset_fma, PolytopeCSC_cacheUpdateCoord_vec}}, "only read!"}},
           }));

      add(new CLIF_Option<PolytopeJIT_Generator>(&PolytopeJIT_generator,'f',"PolytopeJIT_gen","single_rax", {
						     {"single_rax",        {pjit_single_rax,      "single aij at time, load via rax"}},
						     {"single_data",       {pjit_single_data,     "single aij at time, load via data table"}},
						     {"single_data_acc",   {pjit_single_data_acc, "single aij at time, load via data table, more than one acc"}},
						     {"double_data",       {pjit_double_data,     "two aij at time (if possible), load via data table"}},
						     {"quad_data",         {pjit_quad_data,       "four aij at time (if possible), load via data table"}},
						     {"quad_data_acc",     {pjit_quad_data_acc,   "four aij at time (if possible), load via data table, more than one acc"}},
						  }));

      add(new CLIF_Option<intersectCoord_f_t>(&Ellipsoid_T.intersectCoord,'f',"Ellipsoid_intersectCoord","cached_ref", {
                                                     {"ref",        {Ellipsoid_intersectCoord_ref, "no cache (ref)"}},
						     {"cached_ref", {Ellipsoid_intersectCoord_cached_ref, "with cache (ref)"}},
						     {"cached_reord", {Ellipsoid_intersectCoord_cached_reord, "with cache (reord)"}},
						     {"cached_reord2", {Ellipsoid_intersectCoord_cached_reord2, "with cache (reord2)"}},
						     {"cached_reord3", {Ellipsoid_intersectCoord_cached_reord3, "with cache (reord3)"}},
						     {"cached_reord_fma", {Ellipsoid_intersectCoord_cached_reord_fma, "with cache (reord_fma)"}} }));

      add(new CLIF_Option<cacheUpdateCoord_f_t>(&Ellipsoid_T.cacheUpdateCoord,'f',"Ellipsoid_cacheUpdateCoord","ref", {
                                                     {"ref",         {Ellipsoid_cacheUpdateCoord_ref, "cacheUpdateCoord (ref)"}},
						     {"c",           {Ellipsoid_cacheUpdateCoord_c, "cacheUpdateCoord (c)"}},
						     {"fma",           {Ellipsoid_cacheUpdateCoord_fma, "cacheUpdateCoord (fma)"}},
						     {"vec",           {Ellipsoid_cacheUpdateCoord_vec, "cacheUpdateCoord (vec)"}},
						     {"vec_u2",           {Ellipsoid_cacheUpdateCoord_vec_u2, "cacheUpdateCoord (vec_u2)"}},
						     {"vec_u4",           {Ellipsoid_cacheUpdateCoord_vec_u4, "cacheUpdateCoord (vec_u4)"}},
						     {"vec2",           {Ellipsoid_cacheUpdateCoord_vec2, "cacheUpdateCoord (vec2)"}},
						     {"vec2_u2",           {Ellipsoid_cacheUpdateCoord_vec2_u2, "cacheUpdateCoord (vec2_u2)"}},
						     {"vec2_u4",           {Ellipsoid_cacheUpdateCoord_vec2_u4, "cacheUpdateCoord (vec2_u4)"}} }));

      add(new CLIF_Option<cacheReset_f_t>(&Ellipsoid_T.cacheReset,'f',"Ellipsoid_cacheReset","ref", {
                       {"ref",         {Ellipsoid_cacheReset_ref, "cacheReset (ref)"}},
						     {"reord",       {Ellipsoid_cacheReset_reord, "cacheReset (reord)"}},
						     {"fma",       {Ellipsoid_cacheReset_fma, "cacheReset (fma)"}} }));

      add(new CLIF_TrippleOption<rand_init_f_t,rand_f_t,rand256d_f_t>(&rand_init_f, &rand_f, &rand256d_f,'f',"rand_f","std_rand", {
                       {"std_rand",          {{std_init,          {std_rand,         std_rand256d}}, "standard rand"}},
                       {"std_rand_chunked",  {{std_init_chunked,  {std_rand_chunked, std_rand256d}}, "standard rand (chunked)"}},
	               {"sr_rand",           {{sr_init,           {sr_random_uint32, sr_rand256d}},  "shift register rand"}},
                       {"sr_rand_chunked",   {{sr_init_chunked,   {sr_rand_chunked,  sr_rand256d}},  "shift register rand (chunked)"}},
                       {"sr_rand_vec",       {{sr_init_vec,       {sr_rand_vec,      sr_rand256d}},  "shift register rand (vec)"}},
		       {"mt_rand",           {{mt_init,           {mt_rand,          std_rand256d}}, "mersenne twister rand"}} // 256i version missing!
                     }));

      add(new CLIF_Option<rd_0_1_f_t>(&prng_get_random_double_0_1,'f',"rd_0_1","ref", {
                       {"ref",   {prng_get_random_double_0_1, "prng_get_random_double_0_1 (ref)"}},
	               {"fast",  {prng_fast_32_get_random_double_0_1, "prng_get_random_double_0_1 (fast)"}}
                     }));


      // number parameters:
      claimOpt('c',"Algorithm Constants");
      add(new CLIF_OptionNumber<int>(&step_size,'c',"step_size","100000", 100, 1e7));
      add(new CLIF_OptionNumber<int>(&walk_size,'c',"walk_size","1", 1, 1e6));
      add(new CLIF_OptionNumber<int>(&rand_chunk_size,'c',"rand_chunk_size","512", 1, 524288)); // Max: 1 page(4096KB) of doubles
      
      add(new CLIF_Option<int>(&volumeVerbose,'c',"verbose","0", {
                       {"0",  {0, "Nothing"}},
	               {"1",  {1, "important steps"}},
	               {"2",  {2, "also less important steps"}},
	               {"3",  {3, "most things"}},
	               {"4",  {4, "all / heavy debug"}},
                     }));


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
FT read_vinci(string filename, Polytope **P, FT *vol);


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



