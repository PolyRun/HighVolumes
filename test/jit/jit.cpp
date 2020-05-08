#include <iostream>
#include <string>
#include <cassert>

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

   // -------------------------------- end tests

   #ifdef NDEBUG
   std::cout<< "WARNING: DEBUG DISABLED!\n";
   #else
   std::cout<< "TESTS COMPLETE.\n";
   #endif
}
