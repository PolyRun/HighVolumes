#include "intersectCoord.h"

void Pjit_intersectCoord_init_single() {
   double t00 = -FT_MAX;
   double t11 = FT_MAX;
   jit_immediate_via_rax(t00,0);
   jit_immediate_via_rax(t11,1);
}

void Pjit_intersectCoord_body_single(const Polytope* p, const int i, const bool useRax, jit_Table_8** t8) {
   // find relevant entries in column i:
   for(int j=0;j<p->m;j++) {
      FT aij = Polytope_get_a(p,j,i);
      if(aij != 0.0) { // TODO: make epsilon
      
         // d*a = aij
         double aijInv = 1.0/aij;

	 if(useRax) {
            jit_immediate_via_rax(aijInv,4);
	 } else {
	    *t8 = jit_immediate_8_via_data(aijInv,4,*t8);
	 }

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
}

void Pjit_intersectCoord_init_double(jit_Table_16** t16) {
   double t00 = -FT_MAX;
   double t11 = FT_MAX;
   *t16 = jit_immediate_16_via_data(t00,t00,0,*t16);
   *t16 = jit_immediate_16_via_data(t11,t11,1,*t16);
}

void Pjit_intersectCoord_body_double(const Polytope* p, const int i, jit_Table_16** t16) {
   // find all pairs in column i:
   int pairs_max[p->m]; int pairs_max_i = 0;
   int pairs_min[p->m]; int pairs_min_i = 0;

   int state = 0;
   // 0: free
   // 1: last was max pair - want to join?
   // 2: last was min pair - want to join?

   for(int j=0;j<p->m;j++) {
      FT aij = Polytope_get_a(p,j,i);
      if(aij != 0.0) { // TODO: make epsilon
         if(aij < 0.0) {
	    // max
	    if(state==0) {
	       pairs_max[pairs_max_i++] = j;
	       state = 1;
	    } else if(state==1) {
	       state = 0;// join
	    } else {
	       pairs_max[pairs_max_i++] = j;
	       state = 1;
	    }
	 } else {
            // min
	    if(state==0) {
	       pairs_min[pairs_min_i++] = j;
	       state = 2;
	    } else if(state==1) {
	       pairs_min[pairs_min_i++] = j;
	       state = 2;
	    } else {
	       state = 0;// join
	    }
	 }
      } else {
         state = 0;// no pair left in next
      }
   }
   //printf("column %d has %d and %d\n",i,pairs_max_i,pairs_min_i);

   // for now: just a single acc strategy:
   int j = 0;
   while(j < pairs_max_i || j < pairs_min_i) {
      if(j < pairs_max_i) {
         int jj = pairs_max[j];
	 double a0 = Polytope_get_a(p,jj+0,i);
	 double a1 = Polytope_get_a(p,jj+1,i);
	 if(a1 >= 0 || jj+1 >= p->m) {a1 = -1.0/FT_MAX;}// make impotent
	 
	 //printf("max block at: %d %d %f %f\n",i,jj,a0,a1);
         
	 *t16 = jit_immediate_16_via_data(1.0/a0,1.0/a1,4,*t16);
         
	 uint32_t cachej = 8*jj;
	 jit_loadu_xmm(jit_rcx,cachej,3);
         jit_vmulpd_xmm(4,3,2);
	 jit_vmaxpd_xmm(0,2,0);
      }
      if(j < pairs_min_i) {
         int jj = pairs_min[j];
	 double a0 = Polytope_get_a(p,jj+0,i);
	 double a1 = Polytope_get_a(p,jj+1,i);
	 if(a1 <= 0 || jj+1 >= p->m) {a1 = 1.0/FT_MAX;}// make impotent
	 
	 //printf("min block at: %d %d %f %f\n",i,jj,a0,a1);
         
	 *t16 = jit_immediate_16_via_data(1.0/a0,1.0/a1,6,*t16);
         
	 uint32_t cachej = 8*jj;
	 jit_loadu_xmm(jit_rcx,cachej,5);
         jit_vmulpd_xmm(6,5,7);
	 jit_vminpd_xmm(1,7,1);
      }

      j++;
   }
   
   jit_permilpd(0b0101,0,2);
   jit_permilpd(0b0101,1,3);
   jit_vmaxpd_xmm(2,0,0);
   jit_vminpd_xmm(3,1,1);
}

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
   // - have two "accs" for t00,t11
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
   
   // --------------------- set up code facilities:
   jit_Table_8* t8 = NULL;
   jit_Table_16* t16 = NULL;

   // ------------------------------------------- initialize t00,t11
   switch(PolytopeJIT_generator) {
      case pjit_single_rax: case pjit_single_data: {
         Pjit_intersectCoord_init_single();
	 break;
      }
      case pjit_double_data: {
         Pjit_intersectCoord_init_double(&t16);
	 break;
      }
      default: {
	 assert(false && "missing gen code");
	 break;
      }
   }

   // ------------------------------------------ switch case head
   //  assert(p->n < 256 && "if this asserts, then extend for n larger!");
   //  uint8_t nn = p->n;
   //  // 0x83,0xff,xx   cmp    xx,%edi
   //  {const uint8_t instr[] = {0x83,0xff,nn}; jit_push(instr,3);}
   //  // TODO: can we drop this???

   //  // -------------- if bad, jump to end
   //  // 0f 87 xx xx xx xx    	ja  xxxx   -- relative jump
   //  {const uint8_t instr[] = {0x0f,0x87,0,0,0,0}; jit_push(instr,6);}
   //  uint8_t* jump_end = jit_head();// prepare to set L_end here
   
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
      
      switch(PolytopeJIT_generator) {
         case pjit_single_rax: {
            Pjit_intersectCoord_body_single(p,i,true,&t8);
            break;
         }
         case pjit_single_data: {
            Pjit_intersectCoord_body_single(p,i,false,&t8);
            break;
         }
         case pjit_double_data: {
            Pjit_intersectCoord_body_double(p,i,&t16);// TODO implement!
            break;
         }
	 default: {
            assert(false && "missing gen code");
            break;
         }
      }
     
      // jump to end:
      // e9 xx xx xx xx       	jmpq xxxx
      {const uint8_t instr[] = {0xe9,1,1,1,1}; jit_push(instr,5);}
      jump[i] = jit_head();
   }
 
   
   // ---------------------------------------- set L_end:
   //  uint32_t offset = jit_head() - jump_end;
   //  uint64_t offset64 = jit_head() - jump_end;
   //  //printf("jump offset: %d %ld\n",offset, offset64);
   //  jit_write(jump_end-4, (uint8_t*)&offset, 4);
   
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
   
   // -------------------------------- finish up code facilities:
   jit_table_8_consume(t8);
   jit_table_16_consume(t16);

   //jit_print();
}


