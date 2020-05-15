#include "intersectCoord.h"

void PolytopeJIT_generate_intersectCoord_ref(const Polytope *p, PolytopeJIT *o) {
   //jit_print();
   
   o->intersectCoord = (pjit_intersectCoord_f_t)jit_head();
   
   // info about calling convention:
   // d -> %rdi
   // t0 -> %rsi
   // t1 -> %rdx
   // cache -> %rcx

   // required computation:
   // initialize t00, t11
   // jump according to d.
   // calculate intersections
   // jump to end, return values
   
   // ------------------ Performance Analysis:
   // -- mul
   // Port: 1,2
   // lat: 5 (4 on skylake)
   // 2 issued per cycle
   //
   // -- min / max
   // Haswell: 3 lat, 1 issued per cycle
   // Skylake: 4 lat, 3 issued per cycle (or only 2?)
   //  -> add instruction!
   //
   // -- the difficulty:
   // - we have nc = nzA/n entries to look at.
   // - some are negative, some positive -> dependencies on t00,t11 for min,max 
   // - latency: 5+3 for a single entry. +table jump overheads.
   // - best case: execute 2 min/max per 3 cycles.
   // - worst case: 1 per 3 cycles.
   // - single entry per column only: 8 latency + overheads, 2 flops -> less than 0.25 fpc.
   // - filled column: 2 min/max + 2 mul per 3 cycles -> maximum 4/3 fpc.
   // - data: 2 doubles per 3 cycles -> 2*8/3 ~ 5 bytes per cycle
   //
   // Could be improved:
   // - vectorize: but only works if only min/max one vector - require the sorting or rows
   // - 2-packing: fairly easy to get almost 2x speedup
   // - 4-packing: extra overhead to permute over 128 boundary at end. but more flexibility! Depends much more on data.
   //
   // - this does require loading the constants from memory/table at end.
   // - unclear to me is if loading immediate from memory or instruction is faster?

   /// // --------------------- test
   /// //  f2 0f 10 01          	movsd  (%rcx),%xmm0
   /// //  f2 0f 11 02          	movsd  %xmm0,(%rdx)
   /// //  f2 0f 10 41 20       	movsd  0x20(%rcx),%xmm0
   /// //  f2 0f 11 06          	movsd  %xmm0,(%rsi)
   /// {const uint8_t instr[] = {0xf2,0x0f,0x10,0x01}; jit_push(instr,4); }
   /// {const uint8_t instr[] = {0xf2,0x0f,0x11,0x02}; jit_push(instr,4); }
   /// {const uint8_t instr[] = {0xf2,0x0f,0x10,0x41,0x20}; jit_push(instr,5); }
   /// {const uint8_t instr[] = {0xf2,0x0f,0x11,0x06}; jit_push(instr,4); }
   /// { uint8_t instr[] = {0xf3,0xc3}; jit_push(instr,2); }
   /// return;
   
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
   // 4c 63 1c b8          	movslq (%rax,%rdi,4),%r11
   {const uint8_t instr[] = {0x4c,0x63,0x1c,0xb8}; jit_push(instr,4);}
   // 49 01 c3             	add    %rax,%r11
   {const uint8_t instr[] = {0x49,0x01,0xc3}; jit_push(instr,3);}
   // 41 ff e3             	jmpq   *%r11
   {const uint8_t instr[] = {0x41,0xff,0xe3}; jit_push(instr,3);}


   // -------------------------------------------------- Jump table
   jit_allign(4);// allign for array of longs below
   uint8_t* table = jit_head();
   uint32_t L_table = jit_head() - set_table;
   jit_write(set_table-4, (uint8_t*)&L_table,4);
   for(int i=0;i<p->n;i++) {
      {const uint8_t instr[] = {1,1,1,1}; jit_push(instr,4);}
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
	 
	    // d*a = aij
	    //movabs $0xff00ff00ff00ff00,%rax
            {const uint8_t instr[] = {0x48,0xb8}; jit_push(instr,2); }
            double aijInv = 1.0/aij;
	    jit_push((const uint8_t*)&aijInv,8);
            // c4 e1 f9 6e e0       	vmovq  %rax,%xmm4
            {const uint8_t instr[] = {0xc4,0xe1,0xf9,0x6e,0xe0}; jit_push(instr,5);}
	    
	    //printf("A j:%d i:%d a:%f aInv:%f bj:%f\n",j,i,aij, aijInv,bj);

            ///  if(j==4) {
            ///  // debug log out.
            ///  // 0f 28 c4             	movaps %xmm4,%xmm0
            ///  {const uint8_t instr[] = {0x0f,0x28,0xc4}; jit_push(instr,3); }
            ///  // 0f 28 cb             	movaps %xmm3,%xmm1
            ///  {const uint8_t instr[] = {0x0f,0x28,0xcb}; jit_push(instr,3); }
            ///  break;
            ///  }


            ///  // f2 0f 10 a1 00 01 00 	movsd  0x100(%rcx),%xmm4
	    ///  {const uint8_t instr[] = {0xf2,0x0f,0x10,0xa1}; jit_push(instr,4);}
	    ///  uint32_t cachej = 0;//8*j;
	    ///  printf("cachej: %d\n",cachej);
	    ///  jit_push((const uint8_t*)&cachej,4);
            ///  
	    ///  if(j==0) {
            ///  // debug log out.
	    ///  // f2 0f 11 26          	movsd  %xmm4,(%rsi)
            ///  { uint8_t instr[] = {0xf2,0x0f,0x11,0x26}; jit_push(instr,4); }
	    ///  // f2 0f 11 1a          	movsd  %xmm3,(%rdx)
            ///  { uint8_t instr[] = {0xf2,0x0f,0x11,0x1a}; jit_push(instr,4); }
	    ///  
	    ///  { uint8_t instr[] = {0xf3,0xc3}; jit_push(instr,2); }
	    ///  return;
            ///  }

	    // aix = cache[j]   -- %rcx
	    // c5 db 59 91 xx xx xx xx  vmulsd xxxx(%rcx),%xmm4,%xmm2
	    {const uint8_t instr[] = {0xc5,0xdb,0x59,0x91}; jit_push(instr,4);}
	    uint32_t cachej = 8*j;
	    jit_push((const uint8_t*)&cachej,4);

            ///  if(j==4) {
            ///  // debug log out.
            ///  // 0f 28 c2             	movaps %xmm2,%xmm0
            ///  {const uint8_t instr[] = {0x0f,0x28,0xc2}; jit_push(instr,3); }
            ///  // 0f 28 cb             	movaps %xmm3,%xmm1
            ///  {const uint8_t instr[] = {0x0f,0x28,0xcb}; jit_push(instr,3); }
            ///  break;
            ///  }

	    // we already know if min or max!
            if(aij < 0.0) {
	       //printf("do max\n");
               //c5 e9 5f c0          	vmaxpd %xmm0,%xmm2,%xmm0
               {const uint8_t instr[] = {0xc5,0xe9,0x5f,0xc0}; jit_push(instr,4);}
	    } else {
	       //printf("do min\n");
               //c5 e9 5d c9          	vminpd %xmm1,%xmm2,%xmm1
               {const uint8_t instr[] = {0xc5,0xe9,0x5d,0xc9}; jit_push(instr,4);}
	    }
	 }
      }
      
      // jump to end:
      // e9 xx xx xx xx       	jmpq xxxx
      {const uint8_t instr[] = {0xe9,1,1,1,1}; jit_push(instr,5);}
      jump[i] = jit_head();
   }
 
   
   // ---------------------------------------- set L_end:
   uint32_t offset = jit_head() - jump_end;
   uint64_t offset64 = jit_head() - jump_end;
   //printf("jump offset: %d %ld\n",offset, offset64);
   jit_write(jump_end-4, (uint8_t*)&offset, 4);
   
   // make all ends go here
   for(int i=0;i<p->n;i++) {
      uint32_t jump_offset = jit_head()-jump[i];
      uint64_t jump_offset64 = (uint64_t)jit_head()-(uint64_t)jump[i];
      jit_write(jump[i]-4, (uint8_t*)&jump_offset, 4);
      //printf("jump offset: %d %x %lx \n",i,jump_offset,jump_offset64);
      assert(((uint64_t)jump_offset) == jump_offset64);
   }

   // -------------------------------------------- move t00, t11 back
   //f2 0f 11 06          	movsd  %xmm0,(%rsi)
   { uint8_t instr[] = {0xf2,0x0f,0x11,0x06}; jit_push(instr,4); }
   //f2 0f 11 0a          	movsd  %xmm1,(%rdx)
   { uint8_t instr[] = {0xf2,0x0f,0x11,0x0a}; jit_push(instr,4); }
   
   
   // ---- rep ret
   { uint8_t instr[] = {0xf3,0xc3}; jit_push(instr,2); }
   
   o->intersectCoord_bytes = (void*)jit_head() - (void*)o->intersectCoord;
   
   //jit_print();
}


