#include "volume_helper.hpp"
#include "volume_cost.hpp"
#include <cassert>
void xyz_f1_cost(const int n) {
   pc_stack().log(n,2*n, "bogus");
}

void dotProduct_cost_ref(const int n) {
   pc_stack().log(2*n,2*n*sizeof(FT), "dotProduct");
}
void vectorNorm_cost_ref(const int n) {
   pc_stack().log(2*n,n*sizeof(FT), "vectorNorm");
}
void Ball_intersectCoord_cost_ref(const int n) {
   {// frame for vectorNorm
      PC_Frame<vectorNorm_cost_f> frame((void*)vectorNorm);
      frame.costf()(n);
   }

   // read 1
   // div 1
   // add 4
   // mul 9
   // sqrt 1
   pc_stack().log(15,sizeof(FT), "quad. eq.");
}

void volume_cost_ref(const int n, const int bcount, const void** body, const Body_T** type) {
   pc_stack().log(0,n*sizeof(FT), "init_x");
   
   pc_stack().log(0,0, "cacheAlloc");
   pc_stack().log(0,0, "cacheReset - TODO");
   
   pc_stack().log(0,0, "Ball_volume - TODO");
   
   // number of sampling layers (steps)
   size_t l = pc_volume_l; // get from last execution
   size_t s = pc_volume_steps;
   pc_stack().log(0,0, "l layers - TODO");

   {// frame for walk
      PC_Frame<walk_cost_f> frame((void*)walk_f,s);
      frame.costf()(n,bcount,body,type);
   }
}
void walk_cost_ref(const int n, int bcount, const void** body, const Body_T** type) {
   assert(false && "not implemented!");
}

void walkCoord_cost_ref(const int n, int bcount, const void** body, const Body_T** type) {
   int ws = walk_size;
   {// frame for walk_size loop
      PC_Frame_Base loop("loop",ws);

      pc_stack().log(0,0, "random int - TODO");
      
      {// frame for Ball_intersectCoord
         PC_Frame<Ball_intersectCoord_cost_f> frame((void*)Ball_intersectCoord);
         frame.costf()(n);
      }

      pc_stack().log(0,0, "body->intersectCoord - TODO");
      pc_stack().log(0,0, "random double - TODO");
      pc_stack().log(0,0, "body->cacheUpdateCoord - TODO");
   }
}



