#include <iostream>
#include <string>
#include <cassert>

#include <immintrin.h>

extern "C" { // must be included C stlye
#include "../../src/jit/jit.h"
}

int main() {
   std::cout << "start test.\n";

   // -------------------------------- start tests
   printf("memory: %d\n",jit_estimate_npages(1));
   printf("memory: %d\n",jit_estimate_npages(4096 + 1));
   jit_memoryPages();
   jit_memoryPages();
   
   void (*func)();
   func = (void (*)()) jit_head();
    
   std::string msg = "hello I am me!\n";
   size_t msg_size = msg.size();

   #ifdef __linux__
   jit_pushByte(0x48); // mov $1, %rax # set for syscall write
   jit_pushByte(0xc7);
   jit_pushByte(0xc0);
   jit_pushByte(0x01);
   jit_pushByte(0x00);
   jit_pushByte(0x00);
   jit_pushByte(0x00);
   #elif __APPLE__
   jit_pushByte(0x48);
   jit_pushByte(0xc7);
   jit_pushByte(0xc0);
   jit_pushByte(0x04);
   jit_pushByte(0x00);
   jit_pushByte(0x00);
   jit_pushByte(0x02);
   #endif
  
   const uint8_t i_stdout[] = {0x48,0xc7,0xc7,0x01,0x00,0x00,0x00};
   jit_push(i_stdout,7);
   //jit_pushByte(0x48); // mov $1,rdi   # to stdout
   //jit_pushByte(0xc7);
   //jit_pushByte(0xc7);
   //jit_pushByte(0x01);
   //jit_pushByte(0x00);
   //jit_pushByte(0x00);
   //jit_pushByte(0x00);

   jit_pushByte(0x48); // offset to end of code / start of string
   jit_pushByte(0x8d);
   jit_pushByte(0x35);
   jit_pushByte(0x0a); // this is the offset, starts at end of this instr
   jit_pushByte(0x00);
   jit_pushByte(0x00);
   jit_pushByte(0x00);
   
   jit_pushByte(0x48); // mov msg_size, %rdx  # store message size
   jit_pushByte(0xc7);
   jit_pushByte(0xc2);
   jit_pushByte((msg_size & 0xFF) >> 0);// set this
   jit_pushByte((msg_size & 0xFF00) >> 8);// set this
   jit_pushByte((msg_size & 0xFF0000) >> 16);// set this
   jit_pushByte((msg_size & 0xFF000000) >> 24);// set this
   
   jit_pushByte(0x0f); // syscall
   jit_pushByte(0x05);
 
   jit_pushByte(0xc3);  // ret
   
   for(auto c : msg) { // add string data at end
      jit_pushByte(c);
   }

   jit_print();

   func(); // call the function!

   // ------------------------------- simple function with arguments
   {
      int (*func2)(int,int);
      func2 = (int (*)(int,int)) jit_head();
      
      {
         const uint8_t i_mov[] = {0x48,0x89,0xf8}; // mov    %rdi,%rax
         jit_push(i_mov,3);
      }
      {
         const uint8_t i_add[] = {0x48,0x01,0xf0}; // add    %rsi,%rax
         jit_push(i_add,3);
      }
      jit_pushByte(0xc3);  // ret
      
      jit_print();
      
      for(int i=0;i<1000;i++) {
         int a1 = i+5;
         int a2 = i*i;
         int res = func2(a1,a2);
         assert(res==a1+a2);
      }
   }

   // ------------------------------- simple function with double arguments
   {
      double (*func2)(double,double);
      func2 = (double (*)(double,double)) jit_head();
      
      {
         const uint8_t i_add[] = {0xf2,0x0f,0x58,0xc1}; // addsd %xmm1,%xmm0
         jit_push(i_add,4);
      }
      
      jit_pushByte(0xc3);  // ret
      
      jit_print();
      
      for(int i=0;i<1000;i++) {
         double a1 = i+5/2.0;
         double a2 = i*i*1.1;
         double res = func2(a1,a2);
	 std::cout << "i: " << i << " " << a1 << " " << a2 << " " << res << "\n";
         assert(res==a1+a2);
      }
   }

   // ------------------------------- int and double args - operate on conditional
   {
      double (*func2)(int,double);
      func2 = (double (*)(int,double)) jit_head();
      
      {
         const uint8_t i_cmp[] = {0x83,0xff,0x08}; // cmp    $0x8,%edi
         jit_push(i_cmp,3);
      }
      {
         const uint8_t instr[] = {0x77,0x1}; // ja +1  # goes to L_end
         jit_push(instr,2);
      }
      jit_pushByte(0xc3);  // ret
      
      // label: L_end
      {
         const uint8_t instr[] = {0x66,0x0f,0xef,0xc0}; // pxor   %xmm0,%xmm0
         jit_push(instr,4);
      }      
      jit_pushByte(0xc3);  // ret
      
      jit_print();
      
      for(int i=0;i<20;i++) {
         double a = i+5/2.0;
         double res = func2(i,a);
	 std::cout << "cond: " << i << " " << a << " " << res << "\n";
      }
   }

   // ------------------------------- same, just determine relative jump by calculation
   {
      double (*func2)(int,double);
      func2 = (double (*)(int,double)) jit_head();
      
      {
         const uint8_t i_cmp[] = {0x83,0xff,0x08}; // cmp    $0x8,%edi
         jit_push(i_cmp,3);
      }
      {
         const uint8_t instr[] = {0x77,0x0}; // ja +1  # goes to L_end
         jit_push(instr,2);
      }
      uint8_t* jump_end = jit_head();// prepare to set L_end here
      

      jit_pushByte(0xc3);  // ret
      
      *(jump_end-1) = (uint8_t)(jit_head() - jump_end); // got set L_end for jump above
      // label: L_end
      {
         const uint8_t instr[] = {0x66,0x0f,0xef,0xc0}; // pxor   %xmm0,%xmm0
         jit_push(instr,4);
      }      
      jit_pushByte(0xc3);  // ret
      
      jit_print();
       
      for(int i=0;i<20;i++) {
         double a = i+5/2.0;
         double res = func2(i,a);
	 std::cout << "cond: " << i << " " << a << " " << res << "\n";
      }
   }

   // ------------------------------- int -> long table lookup
   {
      long (*func2)(int);
      func2 = (long (*)(int)) jit_head();
      
      {
         const uint8_t i_cmp[] = {0x83,0xff,0x08}; // cmp    $0x8,%edi
         jit_push(i_cmp,3);
      }
      {
         const uint8_t instr[] = {0x77,0x0}; // ja +1  # goes to L_end
         jit_push(instr,2);
      }
      uint8_t* jump_end = jit_head();// prepare to set L_end here

      {
         const uint8_t instr[] = {0x48,0x8d,0x05,0x00,0x00,0x00,0x00}; // lea   L_table(%rip),%rax
         jit_push(instr,7);
      }
      uint8_t* set_table = jit_head();// prepare to set L_table here

      {
         const uint8_t instr[] = {0x89,0xff}; // mov    %edi,%edi
         jit_push(instr,2);
      }
      {
         const uint8_t instr[] = {0x48,0x63,0x04,0xb8}; // movslq (%rax,%rdi,4),%rax
         jit_push(instr,4);
      }
      {
         const uint8_t instr[] = {0xf3,0xc3}; // mov    repz retq
         jit_push(instr,2);
      }
      
      jit_allign(4);// allign for array of longs below

      *(set_table-4) = (uint8_t)(jit_head() - set_table); // got set L_end for jump above
      for(int i=0;i<9;i++) {
         const uint8_t instr[] = {(uint8_t)(i*i+5),0,0,0}; // mov    %edi,%edi
         jit_push(instr,4);
      }

      *(jump_end-1) = (uint8_t)(jit_head() - jump_end); // got set L_end for jump above
      // label: L_end
      {
         const uint8_t instr[] = {0x48,0x31,0xc0}; // xor    %rax,%rax
         jit_push(instr,3);
      }      
      jit_pushByte(0xc3);  // ret
      
      jit_print();
       
      for(int i=-10;i<20;i++) {
         long res = func2(i);
	 std::cout << "table: " << i << " " << res << "\n";
         assert((i<0 && res == 0)  || (i>8 && res==0) || (res==(i*i+5)));
      }
   }

   // ------------------------------- (int,double) -> execute depending on int, return double
   {
      double (*func2)(int,double);
      func2 = (double (*)(int,double)) jit_head();
      
      {
         const uint8_t i_cmp[] = {0x83,0xff,0x08}; // cmp    $0x8,%edi
         jit_push(i_cmp,3);
      }
      {
         const uint8_t instr[] = {0x77,0x0}; // ja +1  # goes to L_end
         jit_push(instr,2);
      }
      uint8_t* jump_end = jit_head();// prepare to set L_end here

      {
         const uint8_t instr[] = {0x48,0x8d,0x15,0x00,0x00,0x00,0x00}; // lea   L_table(%rip),%rdx
         jit_push(instr,7);
      }
      uint8_t* set_table = jit_head();// prepare to set L_table here

      {
         const uint8_t instr[] = {0x89,0xff}; // mov    %edi,%edi
         jit_push(instr,2);
      }
      {
         const uint8_t instr[] = {0x48,0x63,0x04,0xba}; // movslq (%rdx,%rdi,4),%rax
         jit_push(instr,4);
      }
      {
         const uint8_t instr[] = {0x48,0x01,0xd0}; // add    %rdx,%rax
         jit_push(instr,3);
      }
      {
         const uint8_t instr[] = {0xff,0xe0}; // jmpq   *%rax
         jit_push(instr,2);
      }
 
      jit_allign(4);// allign for array of longs below

      uint8_t* table = jit_head();
      *(set_table-4) = (uint8_t)(jit_head() - set_table); // got set L_end for jump above
      for(int i=0;i<9;i++) {
         const uint8_t instr[] = {0,0,0,0}; // mov    %edi,%edi
         jit_push(instr,4);
      }

      for(int i=0;i<9;i++) {// generate code for each, register location in table
         uint8_t* location = jit_head();
	 size_t entry = location - table;
	 std::cout << "entry: " << i << " " << entry << "\n";
	 table[4*i+0] = (entry & 0xff) >> 0;
	 table[4*i+1] = (entry & 0xff00) >> 8;
	 table[4*i+2] = (entry & 0xff0000) >> 16;
	 table[4*i+3] = (entry & 0xff000000) >> 24;
         
	 if(i%3==0) {
	    const uint8_t instr[] = {0x66,0x0f,0xef,0xc0}; // pxor   %xmm0,%xmm0
	    jit_push(instr,4);
	 }
	 jit_pushByte(0xc3);  // ret
      }
      
      size_t offset = jit_head() - jump_end;
      std::cout << "jump offset: " << offset << "\n";
      *(jump_end-1) = (uint8_t)(jit_head() - jump_end); // got set L_end for jump above
      // label: L_end
      {
         const uint8_t instr[] = {0x66,0x0f,0xef,0xc0}; // pxor   %xmm0,%xmm0
         jit_push(instr,4);
      }      
      jit_pushByte(0xc3);  // ret
      
      jit_print();
       
      for(int i=-10;i<20;i++) {
	 double a = i*1.1+2.0;
         double res = func2(i,a);
	 std::cout << "i-call: " << i << " " << a << " " << res << "\n";
      }
   }

   // ------------------------------- get constant
   {
      double (*func2)();
      func2 = (double (*)()) jit_head();
      
      {
         const uint8_t i_instr[] = {0x48,0xbe,0x00,0xff,0x00,0xff,0x00,0xff,0x00,0xff}; // movabs $0xff00ff00ff00ff00,%rsi
         jit_push(i_instr,10);
      }
      uint8_t* p = jit_head();
      double dd = 1.01;
      *(double*)(p-8) = dd;
      {
         const uint8_t i_instr[] = {0xc4,0xe1,0xf9,0x6e,0xc6}; // vmovq  %rsi,%xmm0
         jit_push(i_instr,5);
      }
      jit_pushByte(0xc3);  // ret
      
      jit_print();
       
      double res = func2();
      std::cout << "cal for const: " << res << "\n";
      assert(res == dd);
   }

   {
      std::cout << "CODE GEN EXPERIMENT:\n";
      
      for(int i=0;i<16;i++) {
         jit_clear();
         jit_immediate_via_rax(0.1,i);
	 jit_print();
      }
      std::cout << "immediate_via_rax:\n";
      {
         jit_clear();
         double (*func2)();
         func2 = (double (*)()) jit_head();
         jit_immediate_via_rax(0.1,0);
	 jit_emit_return();
         double res = func2();
         assert(res==0.1);
      }
      std::cout << "immediate_via_data:\n";
      {
         jit_clear();
         double (*func2)();
         func2 = (double (*)()) jit_head();

	 jit_Table_8* t = NULL;// empty list
	 for(int i=0;i<16;i++) {
	    t = jit_immediate_8_via_data(1.0 + i*0.1, i, t);
	 }
	 t = jit_immediate_8_via_data(5.123, 0, t);
	 jit_emit_return();
         jit_table_8_consume(t);
	 double res = func2();
	 assert(res == 5.123);
	 std::cout << "res: " << res << std::endl;
      }
      std::cout << "immediate_16_via_data:\n";
      {
         jit_clear();
         double (*func2)();
         func2 = (double (*)()) jit_head();

	 jit_Table_16* t16 = NULL;// empty list
	 for(int i=0;i<16;i++) {
	    t16 = jit_immediate_16_via_data(4.0,5.0, i, t16);
            jit_permilpd(0b1111, i,i);
	 }
	 jit_emit_return();
         jit_table_16_consume(t16);
	 jit_print();
	 double res = func2();
	 assert(res == 5.0);
	 std::cout << "res: " << res << std::endl;
      }
      std::cout << "permilpd:\n";
      {
         jit_clear();
         double (*func2)();
         func2 = (double (*)()) jit_head();

	 jit_Table_16* t16 = NULL;// empty list
	 t16 = jit_immediate_16_via_data(4.0,5.0, 13, t16);
	 for(int i=12;i>=0;i--) {
            jit_permilpd((i%2)?0b10:0b01, i+1,i);
	 }
	 jit_emit_return();
         jit_table_16_consume(t16);
	 jit_print();
	 double res = func2();
	 assert(res == 5.0);
	 std::cout << "res: " << res << std::endl;
      }
      std::cout << "loadu_xmm - rdi:\n";
      for(int i=0;i<16;i++){
         jit_clear();
         double (*func2)(double*);
         func2 = (double (*)(double*)) jit_head();
         jit_loadu_xmm(jit_rdi,8,i);
         jit_permilpd(0b0101, i,0);
	 jit_emit_return();
	 jit_print();
         double xyz[] = {1.0,2.0,3.0,4.0,5.0};
	 double res = func2((double*)&xyz);
	 assert(res == 3.0);
	 std::cout << "res: " << res << std::endl;
      }
      std::cout << "loadu_xmm - rsi:\n";
      for(int i=0;i<16;i++){
         jit_clear();
         double (*func2)(double*,double*);
         func2 = (double (*)(double*,double*)) jit_head();
         jit_loadu_xmm(jit_rsi,8,i);
         jit_permilpd(0b0101, i,0);
	 jit_emit_return();
	 jit_print();
         double xyz[] = {1.0,2.0,3.0,4.0,5.0};
	 double res = func2(NULL,(double*)&xyz);
	 assert(res == 3.0);
	 std::cout << "res: " << res << std::endl;
      }
      std::cout << "loadu_xmm - rdx:\n";
      for(int i=0;i<16;i++){
         jit_clear();
         double (*func2)(int,int,double*);
         func2 = (double (*)(int,int,double*)) jit_head();
         jit_loadu_xmm(jit_rdx,8,i);
         jit_permilpd(0b0101, i,0);
	 jit_emit_return();
	 jit_print();
         double xyz[] = {1.0,2.0,3.0,4.0,5.0};
	 double res = func2(0,0,(double*)&xyz);
	 assert(res == 3.0);
	 std::cout << "res: " << res << std::endl;
      }
      std::cout << "loadu_xmm - rcx:\n";
      for(int i=0;i<16;i++){
         jit_clear();
         double (*func2)(int,int,int,double*);
         func2 = (double (*)(int,int,int,double*)) jit_head();
         jit_loadu_xmm(jit_rcx,8,i);
         jit_permilpd(0b0101, i,0);
	 jit_emit_return();
	 jit_print();
         double xyz[] = {1.0,2.0,3.0,4.0,5.0};
	 double res = func2(0,0,0,(double*)&xyz);
	 assert(res == 3.0);
	 std::cout << "res: " << res << std::endl;
      }
      
      for(int i=2;i<16;i++) {
         for(int j=2;j<16;j++) {
            for(int k=2;k<16;k++) {
               if(i==j) {continue;}
	       jit_clear();
	       std::cout << "vmulpd xmm " << i << " " << j << " " << k <<"\n";
               double (*func2)(double,double);
               func2 = (double (*)(double,double)) jit_head();
               jit_permilpd(0b1010, 0,i);
               jit_permilpd(0b1010, 1,j);
	       jit_vmulpd_xmm(i,j,k);
               jit_permilpd(0b1010, k,0);
	       jit_emit_return();
	       jit_print();
	       double res = func2(3.0,7.0);
	       assert(res == 21.0);
	    }
	 }
      }
      for(int i=2;i<16;i++) {
         for(int j=2;j<16;j++) {
            for(int k=2;k<16;k++) {
               if(i==j) {continue;}
	       jit_clear();
	       std::cout << "vmaxpd xmm " << i << " " << j << " " << k <<"\n";
               double (*func2)(double,double);
               func2 = (double (*)(double,double)) jit_head();
               jit_permilpd(0b1010, 0,i);
               jit_permilpd(0b1010, 1,j);
	       jit_vmaxpd_xmm(i,j,k);
               jit_permilpd(0b1010, k,0);
	       jit_emit_return();
	       jit_print();
	       double res = func2(3.0,7.0);
	       assert(res == 7.0);
	    }
	 }
      }
      for(int i=2;i<16;i++) {
         for(int j=2;j<16;j++) {
            for(int k=2;k<16;k++) {
               if(i==j) {continue;}
	       jit_clear();
	       std::cout << "vminpd xmm " << i << " " << j << " " << k <<"\n";
               double (*func2)(double,double);
               func2 = (double (*)(double,double)) jit_head();
               jit_permilpd(0b1010, 0,i);
               jit_permilpd(0b1010, 1,j);
	       jit_vminpd_xmm(i,j,k);
               jit_permilpd(0b1010, k,0);
	       jit_emit_return();
	       jit_print();
	       double res = func2(11.0,9.0);
	       assert(res == 9.0);
	    }
	 }
      }

      std::cout << "mul_mem_xmm - rdi:\n";
      for(int i=0;i<16;i++){
         for(int j=0;j<16;j++){
            jit_clear();
            double (*func2)(double*);
            func2 = (double (*)(double*)) jit_head();

	    jit_Table_16* t16 = NULL;// empty list
	    t16 = jit_immediate_16_via_data(4.0,5.0, i, t16);

	    jit_vmulpd_mem_xmm(jit_rdi,8,i,j);
            jit_permilpd(0b0101, j,0);
	    jit_emit_return();
            jit_table_16_consume(t16);
	    jit_print();
            double xyz[] = {1.0,2.0,3.0,4.0,5.0};
	    double res = func2((double*)&xyz);
	    assert(res == 15.0);
	    std::cout << "res: " << res << std::endl;
         }
      }

      std::cout << "mul_mem_xmm - rcx:\n";
      for(int i=0;i<16;i++){
         for(int j=0;j<16;j++){
            jit_clear();
            double (*func2)(int,int,int,double*);
            func2 = (double (*)(int,int,int,double*)) jit_head();

	    jit_Table_16* t16 = NULL;// empty list
	    t16 = jit_immediate_16_via_data(4.0,5.0, i, t16);

	    jit_vmulpd_mem_xmm(jit_rcx,8,i,j);
            jit_permilpd(0b0101, j,0);
	    jit_emit_return();
            jit_table_16_consume(t16);
	    jit_print();
            double xyz[] = {1.0,2.0,3.0,4.0,5.0};
	    double res = func2(0,0,0,(double*)&xyz);
	    assert(res == 15.0);
	    std::cout << "res: " << res << std::endl;
         }
      }

      std::cout << "immediate_32_via_data:\n";
      {
         jit_clear();
         double (*func2)();
         func2 = (double (*)()) jit_head();

	 jit_Table_32* t32 = NULL;// empty list
	 for(int i=0;i<16;i++) {
	    t32 = jit_immediate_32_via_data(4.0,5.0,6.0,7.0, i, t32);
            jit_permpd(0b11111111, i,i);
	 }
	 jit_emit_return();
         jit_table_32_consume(t32);
	 jit_print();
	 double res = func2();
	 assert(res == 7.0);
	 std::cout << "res: " << res << std::endl;
      }

      std::cout << "storeu_xmm:\n";
      {
         jit_clear();
         double (*func2)(double*);
         func2 = (double (*)(double*)) jit_head();
	 int n = 16;
	 double* x = (double*)(aligned_alloc(32, 2*n*sizeof(double))); // align this to 32

	 jit_Table_16* t16 = NULL;// empty list
	 for(int i=0;i<16;i++) {
            int ii = 2*i;
	    t16 = jit_immediate_16_via_data(ii+0,ii+1, i, t16);
	    jit_storeu_xmm(i,jit_rdi,16*i);
	 }
	 jit_emit_return();
         jit_table_16_consume(t16);
	 jit_print();
	 double res = func2(x);
	 for(int i=0;i<16*2;i++) {
            assert(x[i]==i);
	 }
	 std::cout << "res: " << res << std::endl;
      }
 
      std::cout << "storeu_ymm:\n";
      {
         jit_clear();
         double (*func2)(double*);
         func2 = (double (*)(double*)) jit_head();
	 int n = 16;
	 double* x = (double*)(aligned_alloc(32, 4*n*sizeof(double))); // align this to 32

	 jit_Table_32* t32 = NULL;// empty list
	 for(int i=0;i<16;i++) {
            int ii = 4*i;
	    t32 = jit_immediate_32_via_data(ii+0,ii+1,ii+2,ii+3, i, t32);
	    jit_storeu_ymm(i,jit_rdi,32*i);
	 }
	 jit_emit_return();
         jit_table_32_consume(t32);
	 jit_print();
	 double res = func2(x);
	 for(int i=0;i<16*4;i++) {
            assert(x[i]==i);
	 }
	 std::cout << "res: " << res << std::endl;
      }
  
      std::cout << "fma_ymm:\n";
      {
	 double* x = (double*)(aligned_alloc(32, 4*16*sizeof(double))); // align this to 32
         for(int i=0;i<16;i++) {
            for(int j=0;j<16;j++) {
	       if(i==j) {continue;}
	       jit_clear();
               std::cout << "fmadd test ymm " << i << " " << j << "\n";
               double (*func2)(double*);
               func2 = (double (*)(double*)) jit_head();
	       jit_Table_32* t32 = NULL;// empty list
	       
               for(int k=0;k<16;k++) {
		  int ii = k*4;
	          t32 = jit_immediate_32_via_data(ii+0,ii+1,ii+2,ii+3, i, t32);
	          t32 = jit_immediate_32_via_data(3,3,3,3, j, t32);
	          jit_vfmad213pd_mem_ymm(jit_rdi,ii*8,i,j);
		  jit_storeu_ymm(j,jit_rdi,ii*8);
	       }
	       for(int k=0;k<16*4;k++) {x[k]=2.0;}

	       jit_emit_return();
               jit_table_32_consume(t32);
	       double res = func2(x);
	       for(int k=0;k<16*4;k++) {
		  assert(x[k]==2+k*3);
	       }
	    }
	 }
      }
 
      std::cout << "fma_xmm:\n";
      {
	 double* x = (double*)(aligned_alloc(32, 2*16*sizeof(double))); // align this to 32
         for(int i=0;i<16;i++) {
            for(int j=0;j<16;j++) {
	       if(i==j) {continue;}
	       jit_clear();
               std::cout << "fmadd test xmm " << i << " " << j << "\n";
               double (*func2)(double*);
               func2 = (double (*)(double*)) jit_head();
	       jit_Table_16* t16 = NULL;// empty list
	       
               for(int k=0;k<16;k++) {
		  int ii = k*2;
	          t16 = jit_immediate_16_via_data(ii+0,ii+1, i, t16);
	          t16 = jit_immediate_16_via_data(1,1, j, t16);
	          jit_vfmad213pd_mem_xmm(jit_rdi,ii*8,i,j);
		  jit_storeu_xmm(j,jit_rdi,ii*8);
	       }
	       for(int k=0;k<16*2;k++) {x[k]=2.0;}

	       jit_emit_return();
               jit_table_16_consume(t16);
	       double res = func2(x);
	       for(int k=0;k<16*2;k++) {
		  assert(x[k]==2+k);
	       }
	    }
	 }
      }

      std::cout << "jit_vbroadcastsd_ymm:\n";
      {
         for(int i=0;i<16;i++) {
            for(int j=0;j<16;j++) {
	       if(i==j) {continue;}
	       jit_clear();
               std::cout << "jit_vbroadcastsd_ymm " << i << " " << j << "\n";
               double (*func2)();
               func2 = (double (*)()) jit_head();
	       jit_Table_8* t8 = NULL;// empty list
	       
	       t8 = jit_immediate_8_via_data(3, i, t8);
	       jit_vbroadcastsd_ymm(i,j);
               jit_permpd(0b11111111, j,0);
	       jit_emit_return();
               
	       jit_table_8_consume(t8);
	       double res = func2();
	       assert(res == 3);
	    }
	 }
      }

      std::cout << "jit_vbroadcastsd_mem:\n";
      {
         for(int i=0;i<16;i++) {
	       jit_clear();
               std::cout << "jit_vbroadcastsd_mem " << i <<"\n";
               __m256d (*func2)(double*);
               func2 = (__m256d (*)(double*)) jit_head();
	       
	       jit_vbroadcastsd_mem(jit_rdi,8,i);
               jit_permpd(0b11100100, i,0);
	       jit_emit_return();
               
	       double x[2] = {2,3};
	       __m256d res = func2((double*)x);
	       assert(res[0] == 3);
	       assert(res[1] == 3);
	       assert(res[2] == 3);
	       assert(res[3] == 3);
	 }
      }

      std::cout << "jit_broadcast_sd_via_data:\n";
      {
         for(int i=0;i<16;i++) {
	       jit_clear();
               std::cout << "jit_broadcast_sd_via_data " << i << "\n";
               __m256d (*func2)();
               func2 = (__m256d (*)()) jit_head();
	       jit_Table_8* t8 = NULL;// empty list
	       
	       t8 = jit_broadcast_sd_via_data(3,i,t8);

               jit_permpd(0b11100100, i,0);
	       jit_emit_return();
               
	       jit_table_8_consume(t8);
	       __m256d res = func2();
	       assert(res[0] == 3);
	       assert(res[1] == 3);
	       assert(res[2] == 3);
	       assert(res[3] == 3);
	 }
      }

      std::cout << "fnmadd213_ymm:\n";
      {
	 double* x = (double*)(aligned_alloc(32, 4*16*sizeof(double))); // align this to 32
         for(int i=0;i<16;i++) {
            for(int j=0;j<16;j++) {
	       if(i==j) {continue;}
	       jit_clear();
               std::cout << "fnmadd213 test ymm " << i << " " << j << "\n";
               double (*func2)(double*);
               func2 = (double (*)(double*)) jit_head();
	       jit_Table_32* t32 = NULL;// empty list
	       
               for(int k=0;k<16;k++) {
		  int ii = k*4;
	          t32 = jit_immediate_32_via_data(ii+0,ii+1,ii+2,ii+3, i, t32);
	          t32 = jit_immediate_32_via_data(3,3,3,3, j, t32);
	          jit_vfnmad231pd_mem_ymm(jit_rdi,ii*8,i,j);
		  jit_storeu_ymm(j,jit_rdi,ii*8);
	       }
	       for(int k=0;k<16*4;k++) {x[k]=2.0;}

	       jit_emit_return();
               jit_table_32_consume(t32);
	       double res = func2(x);
	       for(int k=0;k<16*4;k++) {
		  assert(x[k]==3-k*2);
	       }
	    }
	 }
      }
 

   }
   // -------------------------------- end tests

   #ifdef NDEBUG
   std::cout<< "WARNING: DEBUG DISABLED!\n";
   #else
   std::cout<< "TESTS COMPLETE.\n";
   #endif
}
