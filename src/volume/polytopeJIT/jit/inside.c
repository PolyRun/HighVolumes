#include "inside.h"

/*
# compile c to asm:
# gcc -S -fverbose-asm -O2 foo.c

# compile asm to machinecode:
# gcc -c hello.s

# dump machinecode to hex: 
# objdump -D hello.o

find test scripts in test/

*/

void PolytopeJIT_generate_inside_ref(const Polytope *p, PolytopeJIT *o) {
   //jit_print();
   
   // info about calling convention:
   // input pointer -> %rdi
   // output boolean -> %eax
   
   // required computation:
   // for each constraint i:
   //    check ai*x <= b

   o->inside = (pjit_inside_f_t)jit_head();
   
   // ----------------- testing:
   // return false;
   // ------------------
   // xorl %eax, %eax
   // ret
   // ------------------
   // 31 c0
   // c3
   //{ uint8_t instr[] = {0x31,0xc0,0xc3}; jit_push(instr,3); }
   
   // ---------------- initialize ret to false
   // xorl %eax, %eax
   { uint8_t instr[] = {0x31,0xc0}; jit_push(instr,2); }

   // ---------------- emmit one check per constraing:
   const int n = p->n;
   const int m = p->m;

   uint8_t* jump[m];
   for(int i=0;i<m;i++) {
      FT* Ai = Polytope_get_Ai(p,i);
      int nonzero = 0;
      for(int j=0;j<n;j++) {
         FT aij = Ai[j];
	 if(aij != 0.0) {// TODO: check if is 1 or -1???
            //printf("Ai %d %d %f\n",i,j,aij);
	    // load immediate to xmm register
            
	    // movabs $0xff00ff00ff00ff00,%rsi
	    {const uint8_t instr[] = {0x48,0xbe}; jit_push(instr,2); }
	    jit_push((const uint8_t*)&aij,8);
            
            uint32_t offset = 8*j;
	    nonzero++;
	    if(nonzero==1) {
	       //c4 e1 f9 6e c6  vmovq  %rsi,%xmm0
	       {const uint8_t instr[] = {0xc4,0xe1,0xf9,0x6e,0xc6}; jit_push(instr,5); }
	       if(offset < (size_t)0x080) {
                  //  f2 0f 59 47 xx  mulsd  xx(%rdi),%xmm0
	          {const uint8_t instr[] = {0xf2,0x0f,0x59,0x47}; jit_push(instr,4); }
	          jit_push((const uint8_t*)&offset,1);
	       } else {
                  //  f2 0f 59 87 xx xx xx xx  mulsd  xx(%rdi),%xmm0
	          {const uint8_t instr[] = {0xf2,0x0f,0x59,0x87}; jit_push(instr,4); }
	          jit_push((const uint8_t*)&offset,4);
	       }
	    } else {
	       //c4 e1 f9 6e ce  vmovq  %rsi,%xmm1
	       {const uint8_t instr[] = {0xc4,0xe1,0xf9,0x6e,0xce}; jit_push(instr,5); }
	       if(offset < (size_t)0x080) {
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
      
      // load b[i]
      FT bi = Polytope_get_b(p,i);
      // movabs $0xff00ff00ff00ff00,%rsi
      {const uint8_t instr[] = {0x48,0xbe}; jit_push(instr,2); }
      jit_push((const uint8_t*)&bi,8);
      //c4 e1 f9 6e ce  vmovq  %rsi,%xmm1
      {const uint8_t instr[] = {0xc4,0xe1,0xf9,0x6e,0xce}; jit_push(instr,5); }
        
      // if( ... < b[i])
      // 66 0f 2e c1 ucomisd %xmm1,%xmm0
      {const uint8_t instr[] = {0x66,0x0f,0x2e,0xc1}; jit_push(instr,4); }
      
      // 0f 87 xx xx xx xx  ja L_end
      {const uint8_t instr[] = {0x0f,0x87,0x0,0x0,0x0,0x0}; jit_push(instr,6); }
      jump[i] = jit_head();
   }
   // ----- return true if one of the conditions was a hit
   //  0f 96 c0  setbe %al
   { uint8_t instr[] = {0x0f,0x96,0xc0}; jit_push(instr,3); }
   
   // go set all the jump addresses:
   uint8_t* L_end = jit_head();
   for(int i=0;i<m;i++) {
      uint32_t offset = L_end - jump[i];
      jit_write(jump[i]-4,(uint8_t*)&offset,4);
   }

   // ---- rep ret
   { uint8_t instr[] = {0xf3,0xc3}; jit_push(instr,2); }
   
   //jit_print();
}

