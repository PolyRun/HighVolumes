#include "intersect.h"

/*
# compile c to asm:
# gcc -S -fverbose-asm -O3 foo.c
# gcc -S -fverbose-asm -march=native -O3 c.c

# compile asm to machinecode:
# gcc -c hello.s

# dump machinecode to hex: 
# objdump -D hello.o

find test scripts in test/

https://www.felixcloutier.com/x86/index.html

*/

void PolytopeJIT_generate_intersect_ref(const Polytope *p, PolytopeJIT *o) {
   jit_print();
   
   // info about calling convention:
   // x -> %rdi
   // d -> %rsi
   // t0 -> %rdx
   // t1 -> %rcx

   // required computation:
   // 2x dotProduct
   // check if d*a is zero/neg/pos -> mask
   // neg: t00 = max(t00,t)
   // pos: t11 = min(t11,t)

   o->intersect = (pjit_intersect_f_t)jit_head();
   
   // --------------------------------- test: t0=x, t1=d
   //f2 0f 10 07   movsd  (%rdi),%xmm0
   //f2 0f 11 02   movsd  %xmm0,(%rdx)
   //f2 0f 10 06   movsd  (%rsi),%xmm0
   //f2 0f 11 01   movsd  %xmm0,(%rcx)
   //{ uint8_t instr[] = {0xf2,0x0f,0x10,0x07}; jit_push(instr,4); }
   //{ uint8_t instr[] = {0xf2,0x0f,0x11,0x02}; jit_push(instr,4); }
   //{ uint8_t instr[] = {0xf2,0x0f,0x10,0x06}; jit_push(instr,4); }
   //{ uint8_t instr[] = {0xf2,0x0f,0x11,0x01}; jit_push(instr,4); }
   

   // -------------- gathering useful assembly ops:
   // vbroadcastsd    %xmm1, %ymm1
   // vfmadd231sd     .LC1(%rip), %xmm4, %xmm2
   //   -- numbers decide which argument goes where and where result goes, variants exist
   // vmulsd  8(%rsi), %xmm2, %xmm2
   // vmulsd  8(%rdi), %xmm0, %xmm0
   // vdivsd  %xmm2, %xmm0, %xmm0
   // vmaxpd  %ymm2, %ymm1, %ymm1
   // vminpd  %ymm0, %ymm1, %ymm0
   // vsubsd  %xmm1, %xmm0, %xmm2
   // vaddsd  %xmm1, %xmm0, %xmm0
   // vzeroupper
   //    -- zeroes out 128 upper bits of all xmm registers, rm dependencies
   // vcmppd  $30, %ymm3, %ymm2, %ymm0
   //    --  predicate, a,b, dst
   // vblendvpd       %ymm1, %ymm2, %ymm3, %ymm4
   //    --  mask, a,b, dst
   
   // ----------------------------- Plan:
   // store t00 in xmm0, t11 in xmm1
   // compute dot product into xmm2 for d*a, xmm3 for x*a
   // compute t into xmm3
   // compute masks into xmm4, xmm5
   // update via min/max t00 and t11
   
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

   
   //c4 e1 f9 6e d0       	vmovq  %rax,%xmm2
   //c4 e1 f9 6e d8       	vmovq  %rax,%xmm3
   //c4 e1 f9 6e e0       	vmovq  %rax,%xmm4
   //c4 e1 f9 6e e8       	vmovq  %rax,%xmm5


   // -------------------------------------------- move t00, t11 back
   //f2 0f 11 02   movsd  %xmm0,(%rdx)
   { uint8_t instr[] = {0xf2,0x0f,0x11,0x02}; jit_push(instr,4); }
   //f2 0f 11 09   movsd  %xmm1,(%rcx)
   { uint8_t instr[] = {0xf2,0x0f,0x11,0x09}; jit_push(instr,4); }

   // ---- rep ret
   { uint8_t instr[] = {0xf3,0xc3}; jit_push(instr,2); }
   
   jit_print();
}

