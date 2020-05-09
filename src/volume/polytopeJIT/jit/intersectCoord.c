#include "intersectCoord.h"

void PolytopeJIT_generate_intersectCoord_ref(const Polytope *p, PolytopeJIT *o) {
   jit_print();
   
   o->intersectCoord = (pjit_intersectCoord_f_t)jit_head();
   { uint8_t instr[] = {0xc3}; jit_push(instr,1); }
   jit_print();
}


