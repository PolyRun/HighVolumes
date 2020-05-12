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
   //jit_print();
   
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
   //   -- a,b, dst    -- dst += a*b
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
   // compute dot product into xmm2 for x*a, xmm3 for d*a
   // compute t into xmm2
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

   
   // -------------------------------------------- for each constraint
   //c4 e1 f9 6e d0       	vmovq  %rax,%xmm2
   //c4 e1 f9 6e d8       	vmovq  %rax,%xmm3
   //c4 e1 f9 6e e0       	vmovq  %rax,%xmm4
   //c4 e1 f9 6e e8       	vmovq  %rax,%xmm5
   for(int i=0;i<p->m;i++) {
      // compute dot product into xmm2 for x*a, xmm3 for d*a
      // c5 e9 57 d2          	vxorpd %xmm2,%xmm2,%xmm2
      {const uint8_t instr[] = {0xc5,0xe9,0x57,0xd2}; jit_push(instr,4);}
      // c5 e1 57 db          	vxorpd %xmm3,%xmm3,%xmm3 
      {const uint8_t instr[] = {0xc5,0xe1,0x57,0xdb}; jit_push(instr,4);}

      FT* Ai = Polytope_get_Ai(p,i);
      for(int j=0;j<p->n;j++) {
         FT aij = Ai[j];
	 if(aij != 0.0) {// TODO: check if is 1 or -1???
            //printf("Ai %d %d %f\n",i,j,aij);
            
	    //movabs $0xff00ff00ff00ff00,%rax
            {const uint8_t instr[] = {0x48,0xb8}; jit_push(instr,2); }
            jit_push((const uint8_t*)&aij,8);
            // c4 e1 f9 6e e0       	vmovq  %rax,%xmm4
            {const uint8_t instr[] = {0xc4,0xe1,0xf9,0x6e,0xe0}; jit_push(instr,5);}
            
            
            // c4 e2 d9 b9 57 08    	vfmadd231sd 0x8(%rdi),%xmm4,%xmm2
            // c4 e2 d9 b9 57 10    	vfmadd231sd 0x10(%rdi),%xmm4,%xmm2
            // c4 e2 d9 b9 97 00 01 	vfmadd231sd 0x100(%rdi),%xmm4,%xmm2
            // 00 00 
            // c4 e2 d9 b9 5e 08    	vfmadd231sd 0x8(%rsi),%xmm4,%xmm3
            // c4 e2 d9 b9 9e 00 01 	vfmadd231sd 0x100(%rsi),%xmm4,%xmm3
            // 00 00

            // vfmadd231sd     a,b,dst
	    // a,b,dst -- dst += a*b
            
            	    
            // c4 e2 d9 b9 97 xx xx xx xx  vfmadd231sd xx(%rdi),%xmm4,%xmm2
	    {const uint8_t instr[] = {0xc4,0xe2,0xd9,0xb9,0x97}; jit_push(instr,5); }
            uint32_t offset = j*8;
	    jit_push((const uint8_t*)&offset,4);
            
	    // c4 e2 d9 b9 9e xx xx xx xx  vfmadd231sd xx(%rsi),%xmm4,%xmm3
	    {const uint8_t instr[] = {0xc4,0xe2,0xd9,0xb9,0x9e}; jit_push(instr,5); }
	    jit_push((const uint8_t*)&offset,4);

	    /// // ----- debug log out.
	    /// // 0f 28 86 00 01 00 00 	movaps 0x100(%rsi),%xmm0
            /// //{const uint8_t instr[] = {0x0f,0x28,0x86,0,0,0,0}; jit_push(instr,5); }
	    /// // 0f 28 06             	movaps (%rsi),%xmm0
	    /// {const uint8_t instr[] = {0x0f,0x28,0x06}; jit_push(instr,3); }
	    /// // 0f 28 cc             	movaps %xmm4,%xmm1
	    /// {const uint8_t instr[] = {0x0f,0x28,0xcc}; jit_push(instr,3); }
	    /// break;
	 }
      }
      /// break;
      
      // if(i==0) {
      // // debug log out.
      // // 0f 28 c2             	movaps %xmm2,%xmm0
      // {const uint8_t instr[] = {0x0f,0x28,0xc2}; jit_push(instr,3); }
      // // 0f 28 cb             	movaps %xmm3,%xmm1
      // {const uint8_t instr[] = {0x0f,0x28,0xcb}; jit_push(instr,3); }
      // break;
      // }

      // assume: xmm2, xmm3 hold x*a, d*a
      
      // load bi into xmm4
      //movabs $0xff00ff00ff00ff00,%rax
      {const uint8_t instr[] = {0x48,0xb8}; jit_push(instr,2); }
      double bi = Polytope_get_b(p,i);
      jit_push((const uint8_t*)&bi,8);
      // c4 e1 f9 6e e0       	vmovq  %rax,%xmm4
      {const uint8_t instr[] = {0xc4,0xe1,0xf9,0x6e,0xe0}; jit_push(instr,5);}
       
      // compute t into xmm2
      // c5 db 5c d2          	vsubsd %xmm2,%xmm4,%xmm2
      {const uint8_t instr[] = {0xc5,0xdb,0x5c,0xd2}; jit_push(instr,4);} 
      // c5 eb 5e d3          	vdivsd %xmm3,%xmm2,%xmm2
      {const uint8_t instr[] = {0xc5,0xeb,0x5e,0xd3}; jit_push(instr,4);} 
      
      /// if(i==2) {
      /// // debug log out.
      /// // 0f 28 c2             	movaps %xmm2,%xmm0
      /// {const uint8_t instr[] = {0x0f,0x28,0xc2}; jit_push(instr,3); }
      /// // 0f 28 cb             	movaps %xmm3,%xmm1
      /// {const uint8_t instr[] = {0x0f,0x28,0xcb}; jit_push(instr,3); }
      /// /// // 0f 28 cc             	movaps %xmm4,%xmm1
      /// /// {const uint8_t instr[] = {0x0f,0x28,0xcc}; jit_push(instr,3); }
      /// break;
      /// }
      
      // compute masks into xmm4, xmm5
      // vcmppd  $30, %ymm3, %ymm2, %ymm0
      //    --  predicate, a,b, dst
      // vblendvpd       %ymm1, %ymm2, %ymm3, %ymm4
      //    --  mask, a,b, dst
      
      // c5 c9 57 f6          	vxorpd %xmm6,%xmm6,%xmm6
      {const uint8_t instr[] = {0xc5,0xc9,0x57,0xf6}; jit_push(instr,4);}

      // vcmppd $17, %xmm6, %xmm3, %xmm4 # 17: OP := _CMP_LT_OQ
      // c5 e3 c2 e6 11       	vcmplt_oqsd %xmm6,%xmm3,%xmm4
      {const uint8_t instr[] = {0xc5,0xe3,0xc2,0xe6,0x11}; jit_push(instr,5);}
      
      // vcmppd $30, %xmm6, %xmm3, %xmm5 # 30: OP := _CMP_GT_OQ
      // c5 e3 c2 ee 1e       	vcmpgt_oqsd %xmm6,%xmm3,%xmm5
      {const uint8_t instr[] = {0xc5,0xe3,0xc2,0xee,0x1e}; jit_push(instr,5);}
      
      /// if(i==2) {
      /// // debug log out.
      /// // 0f 28 c4             	movaps %xmm4,%xmm0
      /// {const uint8_t instr[] = {0x0f,0x28,0xc4}; jit_push(instr,3); }
      /// // 0f 28 cd             	movaps %xmm5,%xmm1
      /// {const uint8_t instr[] = {0x0f,0x28,0xcd}; jit_push(instr,3); }
      /// break;
      /// }
  
      // c4 e3 79 4b e2 40    	vblendvpd %xmm4,%xmm2,%xmm0,%xmm4
      {const uint8_t instr[] = {0xc4,0xe3,0x79,0x4b,0xe2,0x40}; jit_push(instr,6);}
      // c4 e3 71 4b ea 50    	vblendvpd %xmm5,%xmm2,%xmm1,%xmm5
      {const uint8_t instr[] = {0xc4,0xe3,0x71,0x4b,0xea,0x50}; jit_push(instr,6);}

      /// if(i==1) {
      /// // debug log out.
      /// // 0f 28 c4             	movaps %xmm4,%xmm0
      /// {const uint8_t instr[] = {0x0f,0x28,0xc4}; jit_push(instr,3); }
      /// // 0f 28 cd             	movaps %xmm5,%xmm1
      /// {const uint8_t instr[] = {0x0f,0x28,0xcd}; jit_push(instr,3); }
      /// break;
      /// }

      // update via min/max t00 and t11
      // c5 d9 5f c0          	vmaxpd %xmm0,%xmm4,%xmm0
      {const uint8_t instr[] = {0xc5,0xd9,0x5f,0xc0}; jit_push(instr,4);}
      // c5 d1 5d c9          	vminpd %xmm1,%xmm5,%xmm1
      {const uint8_t instr[] = {0xc5,0xd1,0x5d,0xc9}; jit_push(instr,4);}
   }

   // -------------------------------------------- move t00, t11 back
   //f2 0f 11 02   movsd  %xmm0,(%rdx)
   { uint8_t instr[] = {0xf2,0x0f,0x11,0x02}; jit_push(instr,4); }
   //f2 0f 11 09   movsd  %xmm1,(%rcx)
   { uint8_t instr[] = {0xf2,0x0f,0x11,0x09}; jit_push(instr,4); }

   // ---- rep ret
   { uint8_t instr[] = {0xf3,0xc3}; jit_push(instr,2); }
   
   //jit_print();
}

