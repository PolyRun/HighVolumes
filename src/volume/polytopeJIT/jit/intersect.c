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
   // x -> %rdi
   // d -> %rsi
   // t0 -> %rdx
   // t1 -> %rcx

   // required computation:
   //

   o->intersect = (pjit_intersect_f_t)jit_head();
   
   // --------------------------------- test: t0=x, t1=d
   //f2 0f 10 07   movsd  (%rdi),%xmm0
   //f2 0f 11 02   movsd  %xmm0,(%rdx)
   //f2 0f 10 06   movsd  (%rsi),%xmm0
   //f2 0f 11 01   movsd  %xmm0,(%rcx)
   { uint8_t instr[] = {0xf2,0x0f,0x10,0x07}; jit_push(instr,4); }
   { uint8_t instr[] = {0xf2,0x0f,0x11,0x02}; jit_push(instr,4); }
   { uint8_t instr[] = {0xf2,0x0f,0x10,0x06}; jit_push(instr,4); }
   { uint8_t instr[] = {0xf2,0x0f,0x11,0x01}; jit_push(instr,4); }

   // ---- rep ret
   { uint8_t instr[] = {0xf3,0xc3}; jit_push(instr,2); }
   
   jit_print();
}

