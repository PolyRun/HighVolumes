#include "intersect.h"

/*
# compile c to asm:
# gcc -S -fverbose-asm -O2 foo.c

# compile asm to machinecode:
# gcc -c hello.s

# dump machinecode to hex: 
# objdump -D hello.o

find test scripts in test/

*/

void PolytopeJIT_generate_intersect_ref(const Polytope *p, PolytopeJIT *o) {
   jit_print();
   
   // info about calling convention:
   
   // required computation:

   o->intersect = (pjit_inside_f_t)jit_head();
   
   // ---- rep ret
   { uint8_t instr[] = {0xf3,0xc3}; jit_push(instr,1); }
   
   jit_print();
}

