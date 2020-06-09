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


void PolytopeJIT_generate_cacheReset4_ref(const Polytope *p, PolytopeJIT *o) {
   const int n = p->n;
   const int m = p->m;

   //Polytope_print(p);
  
   // info about calling convention:
   // x -> %rdi
   // cache -> %rsi

   o->cacheReset4 = (pjit_cacheReset_f_t)jit_head();
   
   jit_Table_8* t8 = NULL;

   for(int i=0;i<m;i++) {
      FT* Ai = Polytope_get_Ai(p,i);// row

      // read bi into ymm0
      double bi = Polytope_get_b(p,i);
      t8 = jit_broadcast_sd_via_data(bi, 0, t8);
      
      for(int j=0;j<n;j++) {
         FT aij = Ai[j];
         if(aij != 0.0) {// TODO: check if is 1 or -1???
	    t8 = jit_broadcast_sd_via_data(aij, 1, t8);
            
	    // dot -= xj * aij
	    size_t offset = 4*8*j;
            jit_vfnmad231pd_mem_ymm(jit_rdi, offset, 1, 0);
         }
      }
      
      // write %ymm0 to (%rsi)
      uint32_t offset = 8*4*i;
      jit_storeu_ymm(0,jit_rsi,offset); // TODO: use alligned!
   }

   // ---- rep ret
   { uint8_t instr[] = {0xf3,0xc3}; jit_push(instr,2); }

   jit_table_8_consume(t8);
}


void PolytopeJIT_generate_cacheReset8_ref(const Polytope *p, PolytopeJIT *o) {
   const int n = p->n;
   const int m = p->m;

   //Polytope_print(p);
  
   // info about calling convention:
   // x -> %rdi
   // cache -> %rsi

   o->cacheReset8 = (pjit_cacheReset_f_t)jit_head();
   
   jit_Table_8* t8 = NULL;

   for(int i=0;i<m;i++) {
      FT* Ai = Polytope_get_Ai(p,i);// row

      // read bi into ymm0
      double bi = Polytope_get_b(p,i);
      t8 = jit_broadcast_sd_via_data(bi, 0, t8);
      t8 = jit_broadcast_sd_via_data(bi, 1, t8);
      
      for(int j=0;j<n;j++) {
         FT aij = Ai[j];
         if(aij != 0.0) {// TODO: check if is 1 or -1???
	    t8 = jit_broadcast_sd_via_data(aij, 2, t8);
            
	    // dot -= xj * aij
	    size_t offset = 8*8*j;
            jit_vfnmad231pd_mem_ymm(jit_rdi, offset, 2, 0);
            jit_vfnmad231pd_mem_ymm(jit_rdi, offset+4*8, 2, 1);
         }
      }
      
      // write %ymm0 to (%rsi)
      uint32_t offset = 8*8*i;
      jit_storeu_ymm(0,jit_rsi,offset); // TODO: use alligned!
      jit_storeu_ymm(1,jit_rsi,offset+4*8); // TODO: use alligned!
   }

   // ---- rep ret
   { uint8_t instr[] = {0xf3,0xc3}; jit_push(instr,2); }

   jit_table_8_consume(t8);
}



void PolytopeJIT_generate_cacheReset_ref(const Polytope *p, PolytopeJIT *o) {
   //jit_print();
   
   // info about calling convention:
   // x -> %rdi
   // cache -> %rsi

   // required computation:
   // cache = b-A*x

   o->cacheReset = (pjit_cacheReset_f_t)jit_head();
   
   const int n = p->n;
   const int m = p->m;
   for(int i=0;i<m;i++) {
      FT* Ai = Polytope_get_Ai(p,i);

      // read bi into xmm0
      //movabs $0xff00ff00ff00ff00,%rax
      {const uint8_t instr[] = {0x48,0xb8}; jit_push(instr,2); }
      double bi = Polytope_get_b(p,i);
      jit_push((const uint8_t*)&bi,8);
      // c4 e1 f9 6e c0       	vmovq  %rax,%xmm0
      {const uint8_t instr[] = {0xc4,0xe1,0xf9,0x6e,0xc0}; jit_push(instr,5);}

      for(int j=0;j<n;j++) {
         FT aij = Ai[j];
	 if(aij != 0.0) {// TODO: check if is 1 or -1???
            //printf("Ai %d %d %f\n",i,j,aij);
	    // load immediate to xmm register
            
	    // movabs $0xff00ff00ff00ff00,%rax
	    // 48 b8 xxxxxxx 	movabs $0xff00ff00ff00ff00,%rax
	    {const uint8_t instr[] = {0x48,0xb8}; jit_push(instr,2); }
	    jit_push((const uint8_t*)&aij,8);
            
            // TODO: make fmsub out of this:

            size_t offset = 8*j;
	    // c4 e1 f9 6e c8       	vmovq  %rax,%xmm1
	    {const uint8_t instr[] = {0xc4,0xe1,0xf9,0x6e,0xc8}; jit_push(instr,5); }
	    if(offset < (size_t)0x080) {
               //  f2 0f 59 4f xx   mulsd  xx(%rdi),%xmm1
	       {const uint8_t instr[] = {0xf2,0x0f,0x59,0x4f}; jit_push(instr,4); }
	       jit_push((const uint8_t*)&offset,1);
	    } else {
               //  f2 0f 59 8f xx xx xx xx   mulsd  xx(%rdi),%xmm1
	       {const uint8_t instr[] = {0xf2,0x0f,0x59,0x8f}; jit_push(instr,4); }
	       jit_push((const uint8_t*)&offset,4);
	    }

            // f2 0f 5c c1   subsd  %xmm1,%xmm0
	    {const uint8_t instr[] = {0xf2,0x0f,0x5c,0xc1}; jit_push(instr,4); }
	 }
      }
      
      // write %xmm0 to (%rsi)
      // c5 f9 d6 86 xx xx xx xx  vmovq  %xmm0,xx(%rsi)
      {const uint8_t instr[] = {0xc5,0xf9,0xd6,0x86}; jit_push(instr,4); }
      uint32_t offset = 8*i;
      jit_push((const uint8_t*)&offset,4);
   }

   // ---- rep ret
   { uint8_t instr[] = {0xf3,0xc3}; jit_push(instr,2); }
   
   //jit_print();

   // ------------------------------- run functions for 4/8-set
   PolytopeJIT_generate_cacheReset4_ref(p,o);
   PolytopeJIT_generate_cacheReset8_ref(p,o);
}

