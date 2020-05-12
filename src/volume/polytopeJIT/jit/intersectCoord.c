#include "intersectCoord.h"

void PolytopeJIT_generate_intersectCoord_ref(const Polytope *p, PolytopeJIT *o) {
   //jit_print();
   
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


   // ------------------------------------------ switch case head
   assert(p->n < 256 && "if this asserts, then extend for n larger!");
   uint8_t nn = p->n;
   // 0x83,0xff,xx   cmp    xx,%edi
   {const uint8_t instr[] = {0x83,0xff,nn}; jit_push(instr,3);}
   
   // -------------- if bad, jump to end
   // 0f 87 xx xx xx xx    	ja  xxxx   -- relative jump
   {const uint8_t instr[] = {0x0f,0x87,0,0,0,0}; jit_push(instr,6);}
   uint8_t* jump_end = jit_head();// prepare to set L_end here
   
   // ------ get ref to jump table
   // 48 8d 05 xx xx xx xx 	lea    xxxx(%rip),%rax
   {const uint8_t instr[] = {0x48,0x8d,0x05,0,0,0,0}; jit_push(instr,7);}
   uint8_t* set_table = jit_head();// prepare to set L_table here
   // 89 ff                	mov    %edi,%edi
   {const uint8_t instr[] = {0x89,0xff}; jit_push(instr,2);}
   // 48 63 0c b8          	movslq (%rax,%rdi,4),%rcx
   {const uint8_t instr[] = {0x48,0x63,0x0c,0xb8}; jit_push(instr,4);}
   // 48 01 c1             	add    %rax,%rcx
   {const uint8_t instr[] = {0x48,0x01,0xc1}; jit_push(instr,3);}
   // ff e1                	jmpq   *%rcx
   {const uint8_t instr[] = {0xff,0xe1}; jit_push(instr,2);}
   
   // -------------------------------------------------- Jump table
   jit_allign(4);// allign for array of longs below
   uint8_t* table = jit_head();
   uint32_t L_table = jit_head() - set_table;
   jit_write(set_table-4, (uint8_t*)&L_table,4);
   for(int i=0;i<p->n;i++) {
      {const uint8_t instr[] = {0,0,0,0}; jit_push(instr,4);}
   }
   
   // ------------------------------ dump code for each i
   uint8_t* jump[p->n];
   for(int i=0;i<p->n;i++) {// generate code for each, register location in table
      uint8_t* location = jit_head();
      uint32_t entry = location - table;
      jit_write(table+4*i, (uint8_t*)&entry, 4);
      
      // find relevant entries in column:
      for(int j=0;j<p->m;j++) {
         FT aij = Polytope_get_a(p,j,i);
	 if(aij != 0.0) {
	    printf("A %d %d %f\n",j,i,aij);
	 
	    // d*a = aij
	    // read b
            FT bj = Polytope_get_b(p,j);

	    // aix = cache[j]

	    // we already know if min or max!
            if(aij < 0.0) {
	       printf("do max\n");
	    } else {
	       printf("do min\n");
	    }
	 }
      }
      
      // put in meat of the case:
      if(i%2 == 0) {
          // c5 f9 57 c0          	vxorpd %xmm0,%xmm0,%xmm0
	  {const uint8_t instr[] = {0xc5,0xf9,0x57,0xc0}; jit_push(instr,4);}
      } else {
	  // c5 f1 57 c9          	vxorpd %xmm1,%xmm1,%xmm1
	  {const uint8_t instr[] = {0xc5,0xf1,0x57,0xc9}; jit_push(instr,4);}
      }
   
      // jump to end:
      // e9 xx xx xx xx       	jmpq xxxx
      {const uint8_t instr[] = {0xe9,0,0,0,0}; jit_push(instr,5);}
      jump[i] = jit_head();
   }
 
   
   // ---------------------------------------- set L_end:
   uint32_t offset = jit_head() - jump_end;
   printf("jump offset: %d\n",offset);
   jit_write(jump_end-4, (uint8_t*)&offset, 4);
   
   // make all ends go here
   for(int i=0;i<p->n;i++) {
      uint32_t jump_offset = jit_head()-jump[i];
      jit_write(jump[i]-4, (uint8_t*)&jump_offset, 4);
   }

   // -------------------------------------------- move t00, t11 back
   //f2 0f 11 06          	movsd  %xmm0,(%rsi)
   { uint8_t instr[] = {0xf2,0x0f,0x11,0x06}; jit_push(instr,4); }
   //f2 0f 11 0a          	movsd  %xmm1,(%rdx)
   { uint8_t instr[] = {0xf2,0x0f,0x11,0x0a}; jit_push(instr,4); }
   
   
   // ---- rep ret
   { uint8_t instr[] = {0xf3,0xc3}; jit_push(instr,2); }
   
   //jit_print();
}


