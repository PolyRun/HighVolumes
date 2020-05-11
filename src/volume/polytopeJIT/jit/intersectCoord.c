#include "intersectCoord.h"

void PolytopeJIT_generate_intersectCoord_ref(const Polytope *p, PolytopeJIT *o) {
   jit_print();
   
   o->intersectCoord = (pjit_intersectCoord_f_t)jit_head();
   
   // info about calling convention:
   // d -> %rdi
   // t0 -> %rsi
   // t1 -> %rdx

   // required computation:
   // initialize t00, t11
   // jump according to d.
   // calculate intersections
   // jump to end, return values
   
   // ------------------------------------------- initialize t00,t11
   double t00 = -FT_MAX;
   double t11 = FT_MAX;
   //movabs $0xff00ff00ff00ff00,%rax
   {const uint8_t instr[] = {0x48,0xb8}; jit_push(instr,2); }
   jit_push((const uint8_t*)&t00,8);
   // c4 e1 f9 6e c0       	vmovq  %rax,%xmm0
   {const uint8_t instr[] = {0xc4,0xe1,0xf9,0x6e,0xc0}; jit_push(instr,5);}

   //movabs $0xff00ff00ff00ff00,%rax
   {const uint8_t instr[] = {0x48,0xb8}; jit_push(instr,2); }
   jit_push((const uint8_t*)&t11,8);
   // c4 e1 f9 6e c8       	vmovq  %rax,%xmm1
   {const uint8_t instr[] = {0xc4,0xe1,0xf9,0x6e,0xc8}; jit_push(instr,5);}


   
   
   // -------------------------------------------- move t00, t11 back
   //f2 0f 11 06          	movsd  %xmm0,(%rsi)
   { uint8_t instr[] = {0xf2,0x0f,0x11,0x06}; jit_push(instr,4); }
   //f2 0f 11 0a          	movsd  %xmm1,(%rdx)
   { uint8_t instr[] = {0xf2,0x0f,0x11,0x0a}; jit_push(instr,4); }
   
   
   // ---- rep ret
   { uint8_t instr[] = {0xf3,0xc3}; jit_push(instr,2); }
   
   jit_print();
}


