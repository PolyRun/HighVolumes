#include "cacheUpdateCoord.h"

void Pjit_cacheUpdateCoord4_body(const Polytope* p, const int i, jit_Table_8** t8, jit_Table_32** t32) {
   // find relevant entries in column:
   for(int j=0;j<p->m;j++) {
      FT aij = Polytope_get_a(p,j,i);
      if(aij != 0.0) {
	 /// *t8 = jit_immediate_8_via_data(-aij,4,*t8);
         *t8 = jit_broadcast_sd_via_data(-aij, 4, *t8);

         // goal:
         // fmadd: cachej += xmm4 * xmm0
         //
         // now we changed to:
         // cachej -= dx * aid.
         // this we can just make into an addition by taking minus of aid
         // 
         // asm:  xmm4 = cachej + xmm0*xmm4
         // c4 e2 f9 a9 a6 xx xx xx xx  vfmadd213sd 0x100(%rsi),%xmm0,%xmm4
         /// {const uint8_t instr[] = {0xc4,0xe2,0xf9,0xa9,0xa6}; jit_push(instr,5); }
         /// uint32_t cachej = 8*j;
         /// jit_push((const uint8_t*)&cachej,4);
         /// // f2 0f 11 a6 xx xx xx xx     movsd  %xmm4,0x100(%rsi)
         /// {const uint8_t instr[] = {0xf2,0x0f,0x11,0xa6}; jit_push(instr,4); }
         /// jit_push((const uint8_t*)&cachej,4);
         uint32_t cachej = 4*8*j;
	 jit_vfmad213pd_mem_ymm(jit_rsi, cachej, 0, 4);
         jit_storeu_ymm(4,jit_rsi,cachej);
      }
   }

}

void PolytopeJIT_generate_cacheUpdateCoord4_ref(const Polytope *p, PolytopeJIT *o) {
   o->cacheUpdateCoord4 = (pjit_cacheUpdateCoord4_f_t)jit_head();
  
   // info about calling convention:
   // d -> %rdi
   // dx -> %ymm0
   // cache -> %rsi
   
   // --------------------- set up code facilities:
   jit_Table_8* t8 = NULL;
   jit_Table_16* t16 = NULL;
   jit_Table_32* t32 = NULL;

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
      
      Pjit_cacheUpdateCoord4_body(p,i,&t8,&t32);
      
      jit_emit_vzeroupper();// make sure to remove false dependencies!
        
      // ---- rep ret
      { uint8_t instr[] = {0xf3,0xc3}; jit_push(instr,2); }
   }

   // -------------------------------- finish up code facilities:
   jit_table_8_consume(t8);
   jit_table_16_consume(t16);
   jit_table_32_consume(t32);
}



void Pjit_cacheUpdateCoord8_body(const Polytope* p, const int i, jit_Table_8** t8, jit_Table_32** t32) {
   // find relevant entries in column:
   for(int j=0;j<p->m;j++) {
      FT aij = Polytope_get_a(p,j,i);
      if(aij != 0.0) {
	 /// *t8 = jit_immediate_8_via_data(-aij,4,*t8);
         *t8 = jit_broadcast_sd_via_data(-aij, 4, *t8);
         *t8 = jit_broadcast_sd_via_data(-aij, 5, *t8);

         // goal:
         // fmadd: cachej += xmm4 * xmm0
         //
         // now we changed to:
         // cachej -= dx * aid.
         // this we can just make into an addition by taking minus of aid
         // 
         // asm:  xmm4 = cachej + xmm0*xmm4
         // c4 e2 f9 a9 a6 xx xx xx xx  vfmadd213sd 0x100(%rsi),%xmm0,%xmm4
         /// {const uint8_t instr[] = {0xc4,0xe2,0xf9,0xa9,0xa6}; jit_push(instr,5); }
         /// uint32_t cachej = 8*j;
         /// jit_push((const uint8_t*)&cachej,4);
         /// // f2 0f 11 a6 xx xx xx xx     movsd  %xmm4,0x100(%rsi)
         /// {const uint8_t instr[] = {0xf2,0x0f,0x11,0xa6}; jit_push(instr,4); }
         /// jit_push((const uint8_t*)&cachej,4);
         uint32_t cachej = 8*8*j;
	 jit_vfmad213pd_mem_ymm(jit_rsi, cachej, 0, 4);
	 jit_vfmad213pd_mem_ymm(jit_rsi, cachej+4*8, 1, 5);
         jit_storeu_ymm(4,jit_rsi,cachej);
         jit_storeu_ymm(5,jit_rsi,cachej+4*8);
      }
   }
}

void PolytopeJIT_generate_cacheUpdateCoord8_ref(const Polytope *p, PolytopeJIT *o) {
   o->cacheUpdateCoord8 = (pjit_cacheUpdateCoord8_f_t)jit_head();
  
   // info about calling convention:
   // d -> %rdi
   // dx0 -> %ymm0
   // dx1 -> %ymm1
   // cache -> %rsi
   
   // --------------------- set up code facilities:
   jit_Table_8* t8 = NULL;
   jit_Table_16* t16 = NULL;
   jit_Table_32* t32 = NULL;

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
      
      Pjit_cacheUpdateCoord8_body(p,i,&t8,&t32);
      
      jit_emit_vzeroupper();// make sure to remove false dependencies!
        
      // ---- rep ret
      { uint8_t instr[] = {0xf3,0xc3}; jit_push(instr,2); }
   }

   // -------------------------------- finish up code facilities:
   jit_table_8_consume(t8);
   jit_table_16_consume(t16);
   jit_table_32_consume(t32);
}

void Pjit_cacheUpdateCoord_body_single(const Polytope* p, const int i, const bool useRax, jit_Table_8** t8) {
   // find relevant entries in column:
   for(int j=0;j<p->m;j++) {
      FT aij = Polytope_get_a(p,j,i);
      if(aij != 0.0) {
	 if(useRax) {
            jit_immediate_via_rax(-aij,4);
	 } else {
	    *t8 = jit_immediate_8_via_data(-aij,4,*t8);
	 }

         // goal:
         // fmadd: cachej += xmm4 * xmm0
         //
         // now we changed to:
         // cachej -= dx * aid.
         // this we can just make into an addition by taking minus of aid
         // 
         // asm:  xmm4 = cachej + xmm0*xmm4
         // c4 e2 f9 a9 a6 xx xx xx xx  vfmadd213sd 0x100(%rsi),%xmm0,%xmm4
         {const uint8_t instr[] = {0xc4,0xe2,0xf9,0xa9,0xa6}; jit_push(instr,5); }
         uint32_t cachej = 8*j;
         jit_push((const uint8_t*)&cachej,4);
         // f2 0f 11 a6 xx xx xx xx     movsd  %xmm4,0x100(%rsi)
         {const uint8_t instr[] = {0xf2,0x0f,0x11,0xa6}; jit_push(instr,4); }
         jit_push((const uint8_t*)&cachej,4);
      }
   }
}

void Pjit_cacheUpdateCoord_init_quad() {
   // go broadcast xmm to ymm
   jit_vbroadcastsd_ymm(0,0);
}

void Pjit_cacheUpdateCoord_body_quad(const Polytope* p, const int i, jit_Table_32** t32) {
   // find quads:
   int last_quad = -10;
   for(int j=0;j<p->m;j++) {
      FT aij = Polytope_get_a(p,j,i);
      if(aij != 0.0 && last_quad < j-3) {
      	 last_quad = j;
	 double aa[4];
	 for(int k=0;k<4;k++) {
	    if(j+k<p->m) {
	       aa[k] = -Polytope_get_a(p,j+k,i);
	    } else {
	       aa[k] = 0;
	    }
	 }
	 *t32 = jit_immediate_32_via_data(aa[0],aa[1],aa[2],aa[3], 4,*t32);
	 //printf("block at: %d %d: %f %f %f %f\n",i,j,aa[0],aa[1],aa[2],aa[3]);

	 uint32_t cachej = 8*j;
         jit_vfmad213pd_mem_ymm(jit_rsi,cachej,0,4);
	 jit_storeu_ymm(4,jit_rsi,cachej);
      }
   }
   jit_emit_vzeroupper();
}

void Pjit_cacheUpdateCoord_init_double() {
   // go broadcast xmm to xmm
   jit_permilpd_xmm(0b0000,0,0);
}

void Pjit_cacheUpdateCoord_body_double(const Polytope* p, const int i, jit_Table_16** t16) {
   // find quads:
   int last_double = -10;
   for(int j=0;j<p->m;j++) {
      FT aij = Polytope_get_a(p,j,i);
      if(aij != 0.0 && last_double < j-1) {
      	 last_double = j;
	 double aa[2];
	 for(int k=0;k<2;k++) {
	    if(j+k<p->m) {
	       aa[k] = -Polytope_get_a(p,j+k,i);
	    } else {
	       aa[k] = 0;
	    }
	 }
	 *t16 = jit_immediate_16_via_data(aa[0],aa[1], 4,*t16);
	 //printf("block at: %d %d: %f %f %f %f\n",i,j,aa[0],aa[1],aa[2],aa[3]);

	 uint32_t cachej = 8*j;
         jit_vfmad213pd_mem_xmm(jit_rsi,cachej,0,4);
	 jit_storeu_xmm(4,jit_rsi,cachej);
      }
   }
}



void PolytopeJIT_generate_cacheUpdateCoord_ref(const Polytope *p, PolytopeJIT *o) {
   //jit_print();
   
   o->cacheUpdateCoord = (pjit_cacheUpdateCoord_f_t)jit_head();
   
   // info about calling convention:
   // d -> %rdi
   // dx -> %xmm0
   // cache -> %rsi

   // required computation:
   // jump according to d.
   // update relevant cache entries
   // return
   
   // ---------------------- Performance Analysis:
   // -- fmadd:
   // 5 latency (4 on skylake)
   // 2 ports -> 2 issued per cycle
   // 
   // -- movq:
   // Haswell: simd log (3)
   // Skylake: vec alu (3)
   // 
   // -- movsd:
   // Haswell/Skylake: fp mov (1) - bottleneck?
   // - possible dependence if not all 128 bits are written?
   // - probably ok, because writing to mem
   //
   // -- so:
   // - no dependencies
   // - could expect 2 flops per cycle.
   // - this would be 4*8 = 32 bytes/cycle - memory bound!

   ///  // ---------------- test 
   ///  // f2 0f 11 06          	movsd  %xmm0,(%rsi)
   ///  { uint8_t instr[] = {0xf2,0x0f,0x11,0x06}; jit_push(instr,4); }
   ///  // ---- rep ret
   ///  { uint8_t instr[] = {0xf3,0xc3}; jit_push(instr,2); }
   ///  return;

   // --------------------- set up code facilities:
   jit_Table_8* t8 = NULL;
   jit_Table_16* t16 = NULL;
   jit_Table_32* t32 = NULL;

   // ------------------------------------------ switch case head
   //  assert(p->n < 256 && "if this asserts, then extend for n larger!");
   //  uint8_t nn = p->n;
   //  // 0x83,0xff,xx   cmp    xx,%edi
   //  {const uint8_t instr[] = {0x83,0xff,nn}; jit_push(instr,3);}
   //  
   //  // -------------- if bad, jump to end
   //  // 0f 87 xx xx xx xx    	ja  xxxx   -- relative jump
   //  {const uint8_t instr[] = {0x0f,0x87,0,0,0,0}; jit_push(instr,6);}
   //  uint8_t* jump_end = jit_head();// prepare to set L_end here
   
   // -------------------- init: potentially broadcast t:
   switch(PolytopeJIT_generator) {
      case pjit_single_rax:
      case pjit_single_data:
      case pjit_single_data_acc: {
         break;
      }
      case pjit_double_data: {
         Pjit_cacheUpdateCoord_init_double(); // go broadcast xmm to xmm
	 break;
      }
      case pjit_quad_data:
      case pjit_quad_data_acc: {
         Pjit_cacheUpdateCoord_init_quad(); // go broadcast xmm to ymm
         break;
      }
      default: {
         assert(false && "missing gen code");
         break;
      }
   }
    


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
            Pjit_cacheUpdateCoord_body_single(p,i,true,&t8);
            break;
         }
         case pjit_single_data:
         case pjit_single_data_acc: {
            Pjit_cacheUpdateCoord_body_single(p,i,false,&t8);
            break;
         }
         case pjit_double_data: {
            Pjit_cacheUpdateCoord_body_double(p,i,&t16);
            break;
         }
	 case pjit_quad_data:
	 case pjit_quad_data_acc: {
            Pjit_cacheUpdateCoord_body_quad(p,i,&t32);
            break;
         }
	 default: {
            assert(false && "missing gen code");
            break;
         }
      }
        
      // ---- rep ret
      { uint8_t instr[] = {0xf3,0xc3}; jit_push(instr,2); }
   }
   
   //  // ---------------------------------------- set L_end:
   //  uint32_t offset = jit_head() - jump_end;
   //  uint64_t offset64 = jit_head() - jump_end;
   //  //printf("jump offset: %d %ld\n",offset, offset64);
   //  jit_write(jump_end-4, (uint8_t*)&offset, 4);
   //  
   //  // ---- rep ret
   //  { uint8_t instr[] = {0xf3,0xc3}; jit_push(instr,2); }
   
   o->cacheUpdateCoord_bytes = (void*)jit_head() - (void*)o->cacheUpdateCoord;

   // -------------------------------- finish up code facilities:
   jit_table_8_consume(t8);
   jit_table_16_consume(t16);
   jit_table_32_consume(t32);

   //jit_print();

   PolytopeJIT_generate_cacheUpdateCoord4_ref(p, o);
   PolytopeJIT_generate_cacheUpdateCoord8_ref(p, o);
}


