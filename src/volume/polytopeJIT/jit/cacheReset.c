#include "cacheReset.h"

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

void PolytopeJIT_generate_cacheReset_ref(const Polytope *p, PolytopeJIT *o) {
   jit_print();
   
   // info about calling convention:
   // x -> %rdi
   // cache -> %rsi

   // required computation:
   // cache = A*x

   o->cacheReset = (pjit_cacheReset_f_t)jit_head();
   
   const int n = p->n;
   const int m = p->m;
   for(int i=0;i<m;i++) {
      FT* Ai = Polytope_get_Ai(p,i);
      int nonzero = 0;
      for(int j=0;j<n;j++) {
         FT aij = Ai[j];
	 if(aij != 0.0) {// TODO: check if is 1 or -1???
            printf("Ai %d %d %f\n",i,j,aij);
	    // load immediate to xmm register
            
	    // movabs $0xff00ff00ff00ff00,%rax
	    // 48 b8 xxxxxxx 	movabs $0xff00ff00ff00ff00,%rax
	    {const uint8_t instr[] = {0x48,0xb8}; jit_push(instr,2); }
	    jit_push((const uint8_t*)&aij,8);
            
            size_t offset = 8*j;
            nonzero++;
	    if(nonzero==1) {
	       // c4 e1 f9 6e c0       	vmovq  %rax,%xmm0
	       {const uint8_t instr[] = {0xc4,0xe1,0xf9,0x6e,0xc0}; jit_push(instr,5); }
	       if(offset < 0x100) {
                  //  f2 0f 59 47 xx  mulsd  xx(%rdi),%xmm0
	          {const uint8_t instr[] = {0xf2,0x0f,0x59,0x47}; jit_push(instr,4); }
	          jit_push((const uint8_t*)&offset,1);
	       } else {
                  //  f2 0f 59 87 xx xx xx xx  mulsd  xx(%rdi),%xmm0
	          {const uint8_t instr[] = {0xf2,0x0f,0x59,0x87}; jit_push(instr,4); }
	          jit_push((const uint8_t*)&offset,4);
	       }
	    } else {
	       // c4 e1 f9 6e c8       	vmovq  %rax,%xmm1
	       {const uint8_t instr[] = {0xc4,0xe1,0xf9,0x6e,0xc8}; jit_push(instr,5); }
	       if(offset < 0x100) {
                  //  f2 0f 59 4f xx   mulsd  xx(%rdi),%xmm1
	          {const uint8_t instr[] = {0xf2,0x0f,0x59,0x4f}; jit_push(instr,4); }
	          jit_push((const uint8_t*)&offset,1);
	       } else {
                  //  f2 0f 59 8f xx xx xx xx   mulsd  xx(%rdi),%xmm1
	          {const uint8_t instr[] = {0xf2,0x0f,0x59,0x8f}; jit_push(instr,4); }
	          jit_push((const uint8_t*)&offset,4);
	       }

               // f2 0f 58 c1   addsd  %xmm1,%xmm0
	       {const uint8_t instr[] = {0xf2,0x0f,0x58,0xc1}; jit_push(instr,4); }
	    }
	 }
      }
      
      // write %xmm0 to (%rsi)
      // c5 f9 d6 86 xx xx xx xx  vmovq  %xmm0,xx(%rsi)
      {const uint8_t instr[] = {0xc5,0xf9,0xd6,0x86}; jit_push(instr,4); }
      uint32_t offset = 8*i;
      jit_push((const uint8_t*)&offset,4);
   }
 

   //  // -------------------------------------------- move t00, t11 back
   //  //f2 0f 11 02   movsd  %xmm0,(%rdx)
   //  { uint8_t instr[] = {0xf2,0x0f,0x11,0x02}; jit_push(instr,4); }
   //  //f2 0f 11 09   movsd  %xmm1,(%rcx)
   //  { uint8_t instr[] = {0xf2,0x0f,0x11,0x09}; jit_push(instr,4); }

   // ---- rep ret
   { uint8_t instr[] = {0xf3,0xc3}; jit_push(instr,2); }
   
   jit_print();
}

