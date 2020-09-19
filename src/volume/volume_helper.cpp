#include "volume_helper.hpp"


Polytope* Polytope_new_box(int n, FT r) {
   Polytope* p = Polytope_new(n, 2*n);

   for(int i=0; i<n; i++) {// for each dim
      Polytope_set_b(p, i,   r);
      Polytope_set_b(p, i+n, r);
      for(int x=0; x<n; x++) {
         Polytope_set_a(p, i,   x, (x==i)?1:0);
         Polytope_set_a(p, i+n, x, (x==i)?-1:0);
      }
   }

   return p;
}

PolytopeT* PolytopeT_new_box(int n, FT r) {
   PolytopeT* p = PolytopeT_new(n, 2*n);

   for(int i=0; i<n; i++) {// for each dim
      PolytopeT_set_b(p, i,   r);
      PolytopeT_set_b(p, i+n, r);
      for(int x=0; x<n; x++) {
         PolytopeT_set_a(p, i,   x, (x==i)?1:0);
         PolytopeT_set_a(p, i+n, x, (x==i)?-1:0);
      }
   }
   PolytopeT_fix_inv(p);
   return p;
}


void initVolumeFunctions(CLIFunctions &clif) {
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
   
   clif.claimOpt('f',"Algorithm Functions");
 
   clif.add(new CLIF_Option<dotProduct_f_t>(&dotProduct,'f',"dotProduct","2acc", {
                                                  {"ref",  {dotProduct_ref, "ref"}},
                                                  {"2acc", {dotProduct_2acc, "2 accumulators"}},
                                                  {"auto1",{dotProduct_auto1,"auto gen 1"}},
                                                  {"auto2",{dotProduct_auto2,"auto gen 2"}},
     					     {"vec1", {dotProduct_vec1,"vectorize"}} }));
   // no for add_long

   clif.add(new CLIF_Option<squaredNorm_f_t>(&squaredNorm,'f',"squaredNorm","ref", {
     					     {"ref", {squaredNorm_ref, "ref"}} }));
   // no for add_long - there is no choice anyway

   clif.add(new CLIF_TrippleOption<squaredNorm_cached_f_t,squaredNorm_cached_reset_f_t,squaredNorm_cached_update_f_t>(
     	 &squaredNorm_cached,&squaredNorm_cached_reset,&squaredNorm_cached_update,
     	 'f',"squaredNorm_cached","ref", {
                     {"no",        {{squaredNorm_cached_ref,  {squaredNorm_cached_reset_ref, squaredNorm_cached_update_ref }}, "no cacheing (pre-ref)"}},
                     {"ref",       {{squaredNorm_cached_refc, {squaredNorm_cached_reset_refc,squaredNorm_cached_update_refc}}, "cacheing (ref)"}},
     	       	}));
   // no for add_long

   clif.add(new CLIF_DoubleOption<shell_cache_init_f_t, shell_idx_f_t>
       (
        &shell_cache_init, &shell_idx,
        'f', "log_f", "nocache",
        {
         {"nocache", {{shell_cache_init_nocache, shell_idx_nocache}, "no caching (ordinary logs)"}},
         {"ref", {{shell_cache_init_ref, shell_idx_ref}, "linear search in cache"}},
         {"binary", {{shell_cache_init_ref, shell_idx_binary}, "binary search in cache"}},
        }));
   // no for add_long

   clif.add(new CLIF_DoubleOption<walk_f_t,volume_f_t>(&walk_f,&volume,'f',"walk_f","walkCoord_ref", {
                                                  {"walk_ref",          {{walk_ref,               volume_ref},          "random direction walk (ref)"}},
     					     {"walkCoord_ref",     {{walkCoord_ref,          volume_ref},          "coordinate walk (ref)"}},
     					     {"walkCoord_1",       {{walkCoord_coord_single, volume_coord_single}, "coordinate walk, cached squaredNorm"}},
     					     {"walkCoord_4",       {{walkCoord_coord_4,      volume_coord_4},      "coordinate walk, cached squaredNorm, 4-way parallel x"}},
     					     {"walkCoord_8",       {{walkCoord_coord_8,      volume_coord_8},      "coordinate walk, cached squaredNorm, 8-way parallel x"}},
     					     }));
   clif.add_long(new CLIF_DoubleOption<walk_f_t,volume_f_t>(&walk_f,&volume,0,"walk_f","walkCoord_8", {
                                                  {"walk_ref",          {{walk_ref,               volume_ref},          "random direction walk (ref)"}},
     					     {"walkCoord_ref",     {{walkCoord_ref,          volume_ref},          "coordinate walk (ref)"}},
     					     {"walkCoord_1",       {{walkCoord_coord_single, volume_coord_single}, "coordinate walk, cached squaredNorm"}},
     					     {"walkCoord_4",       {{walkCoord_coord_4,      volume_coord_4},      "coordinate walk, cached squaredNorm, 4-way parallel x"}},
     					     {"walkCoord_8",       {{walkCoord_coord_8,      volume_coord_8},      "coordinate walk, cached squaredNorm, 8-way parallel x"}},
     					     }));


   clif.add(new CLIF_Option<intersectCoord_f_t>(&Polytope_T.intersectCoord,'f',"Polytope_intersectCoord","cached_ref", {
                                                  {"ref",        {Polytope_intersectCoord_ref, "no cache (ref)"}},
     					     {"cached_ref", {Polytope_intersectCoord_cached_ref, "with cache (ref)"}} }));
   clif.add_long(new CLIF_Option<intersectCoord_f_t>(&Polytope_T.intersectCoord,0,"Polytope_intersectCoord","cached_ref", {
                                                  {"ref",        {Polytope_intersectCoord_ref, "no cache (ref)"}},
     					     {"cached_ref", {Polytope_intersectCoord_cached_ref, "with cache (ref)"}} }));


   clif.add(new CLIF_TrippleOption<intersectCoord_f_t,cacheReset_f_t,cacheUpdateCoord_f_t>(
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
   clif.add_long(new CLIF_TrippleOption<intersectCoord_f_t,cacheReset_f_t,cacheUpdateCoord_f_t>(
     	 &PolytopeT_T.intersectCoord,&PolytopeT_T.cacheReset,&PolytopeT_T.cacheUpdateCoord,
     	 0,"PolytopeT_intersectCoord","cached_b_inv_vec", {
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


   clif.add(new CLIF_TrippleOption<intersectCoord_f_t, cacheReset_f_t, cacheUpdateCoord_f_t>
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
   clif.add_long(new CLIF_TrippleOption<intersectCoord_f_t, cacheReset_f_t, cacheUpdateCoord_f_t>
       (&PolytopeCSC_T.intersectCoord, &PolytopeCSC_T.cacheReset, &PolytopeCSC_T.cacheUpdateCoord,
        0, "PolytopeCSC_intersectCoord", "cached_b_vec_nan_inv",
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


   clif.add(new CLIF_Option<PolytopeJIT_Generator>(&PolytopeJIT_generator,'f',"PolytopeJIT_gen","single_rax", {
     					     {"single_rax",        {pjit_single_rax,      "single aij at time, load via rax"}},
     					     {"single_data",       {pjit_single_data,     "single aij at time, load via data table"}},
     					     {"single_data_acc",   {pjit_single_data_acc, "single aij at time, load via data table, more than one acc"}},
     					     {"double_data",       {pjit_double_data,     "two aij at time (if possible), load via data table"}},
     					     {"quad_data",         {pjit_quad_data,       "four aij at time (if possible), load via data table"}},
     					     {"quad_data_acc",     {pjit_quad_data_acc,   "four aij at time (if possible), load via data table, more than one acc"}},
     					  }));
   clif.add_long(new CLIF_Option<PolytopeJIT_Generator>(&PolytopeJIT_generator,0,"PolytopeJIT_gen","single_data", {
     					     {"single_rax",        {pjit_single_rax,      "single aij at time, load via rax"}},
     					     {"single_data",       {pjit_single_data,     "single aij at time, load via data table"}},
     					     {"single_data_acc",   {pjit_single_data_acc, "single aij at time, load via data table, more than one acc"}},
     					     {"double_data",       {pjit_double_data,     "two aij at time (if possible), load via data table"}},
     					     {"quad_data",         {pjit_quad_data,       "four aij at time (if possible), load via data table"}},
     					     {"quad_data_acc",     {pjit_quad_data_acc,   "four aij at time (if possible), load via data table, more than one acc"}},
     					  }));


   clif.add(new CLIF_Option<intersectCoord_f_t>(&Ellipsoid_T.intersectCoord,'f',"Ellipsoid_intersectCoord","cached_ref", {
                                                  {"ref",        {Ellipsoid_intersectCoord_ref, "no cache (ref)"}},
     					     {"cached_ref", {Ellipsoid_intersectCoord_cached_ref, "with cache (ref)"}},
     					     {"cached_reord", {Ellipsoid_intersectCoord_cached_reord, "with cache (reord)"}},
     					     {"cached_reord2", {Ellipsoid_intersectCoord_cached_reord2, "with cache (reord2)"}},
     					     {"cached_reord3", {Ellipsoid_intersectCoord_cached_reord3, "with cache (reord3)"}},
     					     {"cached_reord_fma", {Ellipsoid_intersectCoord_cached_reord_fma, "with cache (reord_fma)"}} }));
   clif.add_long(new CLIF_Option<intersectCoord_f_t>(&Ellipsoid_T.intersectCoord,0,"Ellipsoid_intersectCoord","cached_reord_fma", {
                                                  {"ref",        {Ellipsoid_intersectCoord_ref, "no cache (ref)"}},
     					     {"cached_ref", {Ellipsoid_intersectCoord_cached_ref, "with cache (ref)"}},
     					     {"cached_reord", {Ellipsoid_intersectCoord_cached_reord, "with cache (reord)"}},
     					     {"cached_reord2", {Ellipsoid_intersectCoord_cached_reord2, "with cache (reord2)"}},
     					     {"cached_reord3", {Ellipsoid_intersectCoord_cached_reord3, "with cache (reord3)"}},
     					     {"cached_reord_fma", {Ellipsoid_intersectCoord_cached_reord_fma, "with cache (reord_fma)"}} }));


   clif.add(new CLIF_Option<cacheUpdateCoord_f_t>(&Ellipsoid_T.cacheUpdateCoord,'f',"Ellipsoid_cacheUpdateCoord","ref", {
                                                  {"ref",         {Ellipsoid_cacheUpdateCoord_ref, "cacheUpdateCoord (ref)"}},
     					     {"c",           {Ellipsoid_cacheUpdateCoord_c, "cacheUpdateCoord (c)"}},
     					     {"fma",           {Ellipsoid_cacheUpdateCoord_fma, "cacheUpdateCoord (fma)"}},
     					     {"vec",           {Ellipsoid_cacheUpdateCoord_vec, "cacheUpdateCoord (vec)"}},
     					     {"vec_u2",           {Ellipsoid_cacheUpdateCoord_vec_u2, "cacheUpdateCoord (vec_u2)"}},
     					     {"vec_u4",           {Ellipsoid_cacheUpdateCoord_vec_u4, "cacheUpdateCoord (vec_u4)"}},
     					     {"vec2",           {Ellipsoid_cacheUpdateCoord_vec2, "cacheUpdateCoord (vec2)"}},
     					     {"vec2_u2",           {Ellipsoid_cacheUpdateCoord_vec2_u2, "cacheUpdateCoord (vec2_u2)"}},
     					     {"vec2_u4",           {Ellipsoid_cacheUpdateCoord_vec2_u4, "cacheUpdateCoord (vec2_u4)"}} }));
   clif.add_long(new CLIF_Option<cacheUpdateCoord_f_t>(&Ellipsoid_T.cacheUpdateCoord,0,"Ellipsoid_cacheUpdateCoord","vec2_u2", {
                                                  {"ref",         {Ellipsoid_cacheUpdateCoord_ref, "cacheUpdateCoord (ref)"}},
     					     {"c",           {Ellipsoid_cacheUpdateCoord_c, "cacheUpdateCoord (c)"}},
     					     {"fma",           {Ellipsoid_cacheUpdateCoord_fma, "cacheUpdateCoord (fma)"}},
     					     {"vec",           {Ellipsoid_cacheUpdateCoord_vec, "cacheUpdateCoord (vec)"}},
     					     {"vec_u2",           {Ellipsoid_cacheUpdateCoord_vec_u2, "cacheUpdateCoord (vec_u2)"}},
     					     {"vec_u4",           {Ellipsoid_cacheUpdateCoord_vec_u4, "cacheUpdateCoord (vec_u4)"}},
     					     {"vec2",           {Ellipsoid_cacheUpdateCoord_vec2, "cacheUpdateCoord (vec2)"}},
     					     {"vec2_u2",           {Ellipsoid_cacheUpdateCoord_vec2_u2, "cacheUpdateCoord (vec2_u2)"}},
     					     {"vec2_u4",           {Ellipsoid_cacheUpdateCoord_vec2_u4, "cacheUpdateCoord (vec2_u4)"}} }));


   clif.add(new CLIF_Option<cacheReset_f_t>(&Ellipsoid_T.cacheReset,'f',"Ellipsoid_cacheReset","ref", {
                    {"ref",         {Ellipsoid_cacheReset_ref, "cacheReset (ref)"}},
     					     {"reord",       {Ellipsoid_cacheReset_reord, "cacheReset (reord)"}},
     					     {"fma",       {Ellipsoid_cacheReset_fma, "cacheReset (fma)"}} }));
   clif.add_long(new CLIF_Option<cacheReset_f_t>(&Ellipsoid_T.cacheReset,0,"Ellipsoid_cacheReset","fma", {
                    {"ref",         {Ellipsoid_cacheReset_ref, "cacheReset (ref)"}},
     					     {"reord",       {Ellipsoid_cacheReset_reord, "cacheReset (reord)"}},
     					     {"fma",       {Ellipsoid_cacheReset_fma, "cacheReset (fma)"}} }));


   clif.add(new CLIF_TrippleOption<rand_init_f_t,rand_f_t,rand256d_f_t>(&rand_init_f, &rand_f, &rand256d_f,'f',"rand_f","std_rand", {
                    {"std_rand",          {{std_init,          {std_rand,         std_rand256d}}, "standard rand"}},
                    {"std_rand_chunked",  {{std_init_chunked,  {std_rand_chunked, std_rand256d}}, "standard rand (chunked)"}},
                    {"sr_rand",           {{sr_init,           {sr_random_uint32, sr_rand256d}},  "shift register rand"}},
                    {"sr_rand_chunked",   {{sr_init_chunked,   {sr_rand_chunked,  sr_rand256d}},  "shift register rand (chunked)"}},
                    {"sr_rand_vec",       {{sr_init_vec,       {sr_rand_vec,      sr_rand256d}},  "shift register rand (vec)"}},
     	       {"mt_rand",           {{mt_init,           {mt_rand,          std_rand256d}}, "mersenne twister rand"}} // 256i version missing!
                  }));
   clif.add_long(new CLIF_TrippleOption<rand_init_f_t,rand_f_t,rand256d_f_t>(&rand_init_f, &rand_f, &rand256d_f,0,"rand_f","sr_rand_vec", {
                    {"std_rand",          {{std_init,          {std_rand,         std_rand256d}}, "standard rand"}},
                    {"std_rand_chunked",  {{std_init_chunked,  {std_rand_chunked, std_rand256d}}, "standard rand (chunked)"}},
                    {"sr_rand",           {{sr_init,           {sr_random_uint32, sr_rand256d}},  "shift register rand"}},
                    {"sr_rand_chunked",   {{sr_init_chunked,   {sr_rand_chunked,  sr_rand256d}},  "shift register rand (chunked)"}},
                    {"sr_rand_vec",       {{sr_init_vec,       {sr_rand_vec,      sr_rand256d}},  "shift register rand (vec)"}},
     	       {"mt_rand",           {{mt_init,           {mt_rand,          std_rand256d}}, "mersenne twister rand"}} // 256i version missing!
                  }));

   clif.add(new CLIF_Option<rd_0_1_f_t>(&prng_get_random_double_0_1,'f',"rd_0_1","ref", {
                    {"ref",   {prng_get_random_double_0_1, "prng_get_random_double_0_1 (ref)"}},
                    {"fast",  {prng_fast_32_get_random_double_0_1, "prng_get_random_double_0_1 (fast)"}}
                  }));
   clif.add_long(new CLIF_Option<rd_0_1_f_t>(&prng_get_random_double_0_1,0,"rd_0_1","fast", {
                    {"ref",   {prng_get_random_double_0_1, "prng_get_random_double_0_1 (ref)"}},
                    {"fast",  {prng_fast_32_get_random_double_0_1, "prng_get_random_double_0_1 (fast)"}}
                  }));


   // number parameters:
   clif.claimOpt('c',"Algorithm Constants");
   clif.add(new CLIF_OptionNumber<int>(&step_size,'c',"step_size","100000", 100, 1e7));
   clif.add_long(new CLIF_OptionNumber<int>(&step_size,0,"step_size","100000", 100, 1e7,
			   "Number of steps per zone. Set higher to increase precision and runtime."));
   
   clif.add(new CLIF_OptionNumber<int>(&walk_size,'c',"walk_size","1", 1, 1e6));
   clif.add_long(new CLIF_OptionNumber<int>(&walk_size,0,"walk_size","1", 1, 1e6, "Number of walks per step (sample point). Not much benefit if >1."));
   
   clif.add(new CLIF_OptionNumber<int>(&rand_chunk_size,'c',"rand_chunk_size","512", 1, 524288)); // Max: 1 page(4096KB) of doubles
   // no for add_long

   clif.add(new CLIF_Option<int>(&volumeVerbose,'c',"verbose","0", {
                    {"0",  {0, "Nothing"}},
                    {"1",  {1, "important steps"}},
                    {"2",  {2, "also less important steps"}},
                    {"3",  {3, "most things"}},
                    {"4",  {4, "all / heavy debug"}},
                  }));
   clif.add_long(new CLIF_Option<int>(&volumeVerbose,'v',"verbose","0", {
                    {"0",  {0, "Nothing"}},
                    {"1",  {1, "important steps"}},
                    {"2",  {2, "also less important steps"}},
                    {"3",  {3, "most things"}},
                    {"4",  {4, "all / heavy debug"}},
                  }));
}


void make_random_poly(const std::vector<double> &ell, int m, Polytope **ret){

  // create n random hyperplanes around ellipsoid ell
  int n = ell.size();
  *ret = Polytope_new(n, m);
  
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  std::normal_distribution<double> distribution(0,1);

  
  for (int i = 0; i < m; i++){

    // sample uar on unit dim-sphere
    std::vector<double> pnt(n);
    double rad = 0;

    for (int j = 0; j < n; j++) {
      double number = distribution(generator);
      rad += number * number;
      pnt[j] = number;
    }
    rad = sqrt(rad);

    
    // normalize point (on surface) and map point to ellipsoid
    // (note it's no longer uniform at that point but maybe good enough, we could make it uniform by discarding each point with a probability proportional to the distortion factor of the matrix inducing the ellipsoid at that point)
    for (int j = 0; j < n; j++) {
        pnt[j] *= ell[j]/rad;
    }

    // add the tangent plane of ellispoid at pnt to poly
    // note that gradient of the ellipsoid at x is [2*ell_i x_i]_{i = 1 to dim}
    double sum = 0;
    for (int j = 0; j < n; j++){
        Polytope_set_a(*ret, i, j, 2*pnt[j] / (ell[j] * ell[j]));
        sum += Polytope_get_a(*ret, i, j) * pnt[j];
    }
    Polytope_set_b(*ret, i, sum);    
  }

}



void polyvest_convert(Polytope *P, vol::Polyvest_p *Q){

    int n = P->n;
    int m = P->m;  

    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            Q->matA(Polytope_get_a(P, i, j), i, j);
        }
        Q->vecb(Polytope_get_b(P, i), i);
    }
  
}


int read_polyvest_p(string filename, Polytope **P){

    ifstream file;
    file.open(filename);


    if (!file.is_open()){
        printf("failed to read polytope");
        return 1;
    }

    int n, m;
    file >> m >> n;

    *P = Polytope_new(n, m);

    FT num;
    for (int i = 0; i < m; i++){
        file >> num;
        Polytope_set_b(*P, i, num);
        for (int j = 0; j < n; j++){
            file >> num;
            Polytope_set_a(*P, i, j, num);
        }
    }

    return 0;
}


// copied from vinci
FT sread_rational_value (char *s);
FT sread_rational_value (char *s) {
    
   char *numerator_s, *denominator_s = NULL, *position, token;
   int sign = 1, i;
   FT numerator, denominator;

   /* determine the sign of the number */
   numerator_s = s;
   if (s [0] == '-')
   {  sign = -1;
      numerator_s++;
   }
   else if (s [0] == '+')
      numerator_s++;

   /* look for a sign '/' and in this case split the number in numerator and denominator */
   position = strchr (numerator_s, '/');
   if (position != NULL)
   {  *position = '\0'; /* terminates the numerator */
      denominator_s = position + 1;
   };

   /* determine the floating point values of numerator and denominator */
   numerator = 0;
   for (i = 0; i < strlen (numerator_s); i++)
   {  token = numerator_s [i];
      if (strchr ("0123456789", token)) /* token is a cypher */
         numerator = 10 * numerator + (int) token - 48;
   }

   if (position != NULL)
   {  denominator = 0;
      for (i = 0; i < strlen (denominator_s); i++)
      {  token = denominator_s [i];
         if (strchr ("0123456789", token)) /* token is a cypher */
            denominator = 10 * denominator + (int) token - 48;
      }
   }
   else denominator = 1;

   return sign * numerator / denominator;
}


FT read_vinci_nr(string in, string type){

    istringstream instr(in);
    if (!type.compare("integer")){
        FT d;
        instr >> d;
        return d;
    }
    else if (!type.compare("real")){
        FT f;
        instr >> f;
        return f;
    }
    else if (!type.compare("rational")){
        
        std::regex rgx("(.*)/(.*)");
        std::smatch matches;

        if(std::regex_search(in, matches, rgx)) {

            FT num, den;
            istringstream ns(matches[1]);
            istringstream ds(matches[2]);
            ns >> num;
            ds >> den;
            //cout << "print " << num/den << " for " << in << "\n";
            return num/den;
        }
        else {
            FT d;
            istringstream s(in);
            s >> d;
            return d;
        }
    }

    assert(0 && "not a valid number");
    return 0.0;

    

}


int read_vinci(string filename, Polytope **P, FT *vol){
    ifstream file;
    file.open(filename);

    if (!file.is_open()){
        printf("failed to read polytope");
        return 1;
    }

    std::string line, type, in;
    std::string b = "begin";
    bool found = false;
    while (std::getline(file, line)){
        if (!b.compare(line)) {
            found = true;
            break;
        }
    }

    if (!found){
        printf("no begin found in vinci file\n");
        return 1;
    }


    int m, n;
    file >> m >> n >> type;
    n--;

    *P = Polytope_new(n,m);
                      
    for (int i = 0; i < m; i++){
        file >> in;
        Polytope_set_b(*P, i, read_vinci_nr(in, type));
        for (int j = 0; j < n; j++){
            file >> in;
            Polytope_set_a(*P, i, j, read_vinci_nr(in, type));            
        }
    }

    // extract solved value
    while (std::getline(file, line)){
        std::regex rgx("Volume:(.*)");
        std::smatch matches;

        if(std::regex_search(line, matches, rgx)) {

            istringstream os(matches[0].str());
            string s;
            os >> s >> *vol;
            return 0;
        }
    }

    *vol = 0;
    return 0;
    
}

int polytope_count_pairs(Polytope* p) {
    const int n= p->n;
    const int m= p->m;

    int count = 0;
    int nz = 0;
    for(int i=0;i<n;i++) {
       double last = 0;
       for(int j=0;j<m;j++) {
          double aij = Polytope_get_a(p,j,i);
	  if(aij!=0) {
             nz++;
	     if(last!=0 && aij*last > 0) {count++;}
	  }
	  if(last!=0) {last =0;} else {last = aij;}
       }
    }
    std::cout << "NZ: " << nz << ", pairs: "<< count << "\n";
    return count;
}

struct OElement {
   int u,v;
   int overlap;
};

Polytope* optimize_polytope(Polytope *p) {
   polytope_count_pairs(p);
   Polytope* q = (Polytope*)Polytope_T.clone(p);
   
   // generator element for each pair of constraints:
   const int n= p->n;
   const int m= p->m;
   std::vector<OElement> oel;
   for(int u=0;u<m;u++) {
      for(int v=0;v<u;v++) {
         int overlap = 0;
	 for(int i=0;i<n;i++) {
	    FT ui = Polytope_get_a(p,u,i);
	    FT vi = Polytope_get_a(p,v,i);
	    if(ui*vi>0) {overlap++;}
	 }
	 if(overlap > 0) {oel.push_back({u,v,overlap});}
      }
   }
   std::sort(oel.begin(),oel.end(), [](const OElement &a,const OElement &b){return a.overlap > b.overlap;});
   //std::cout << "oel " << oel.size() <<  "\n";
   
   std::vector<int> neighbor_1(m,-1);
   std::vector<int> neighbor_2(m,-1);
   evp::Union_find uf(m);
   for(const OElement &el : oel) {
      // check if even possible:
      if(neighbor_2[el.u]!=-1 || neighbor_2[el.v]!=-1) {continue;}
      // check for cycles:
      if(uf.find(el.u) == uf.find(el.v)) {continue;}
      uf.Union(el.u,el.v);

      // allocate u:
      if(neighbor_1[el.u]==-1) {
         neighbor_1[el.u] = el.v;
      } else {
         neighbor_2[el.u] = el.v;
      }
      // allocate v:
      if(neighbor_1[el.v]==-1) {
         neighbor_1[el.v] = el.u;
      } else {
         neighbor_2[el.v] = el.u;
      }
   }
   
   ///  // print check:
   ///  for(int i=0;i<m;i++) {
   ///     std::cout << "n: " << i << " " << neighbor_1[i] << " " << neighbor_2[i] << "\n";
   ///  }

   // go traverse:
   std::vector<bool> marked(m,false);
   std::vector<int> permutation;
   for(int i=0;i<m;i++) {
      //std::cout << "n: " << i << " " << neighbor_1[i] << " " << neighbor_2[i] << "\n";
      if(!marked[i] && neighbor_2[i]==-1) {
         // start here!
         int current = i;
	 int last = -1;
	 while(true) {
            permutation.push_back(current);
	    assert(!marked[current]);
	    marked[current] = true;
	    int next = (neighbor_2[current]==last)?neighbor_1[current]:neighbor_2[current];
	    //std::cout << "push " << current << ", next: "<< next<< "\n";
	    if(next==-1) {break;}
	    last = current;
	    current = next;
	 }
      }
   }

   ///  // print check:
   ///  for(int i=0;i<m;i++) {
   ///     std::cout << "fin: " << i << " " << marked[i] << " " << neighbor_1[i] << " " << neighbor_2[i] << "\n";
   ///  }
   //std::cout << "perm " << permutation.size() << " " << m << "\n";
   assert(permutation.size()==m);
   
   // copy rows according to permutation:
   for(int i=0;i<n;i++) {
      double last = 0;
      for(int j=0;j<m;j++) { 
	 double aij = Polytope_get_a(p,permutation[j],i);
         Polytope_set_a(q,j,i,aij);
      }
   }

   polytope_count_pairs(q);
   return q;
}




