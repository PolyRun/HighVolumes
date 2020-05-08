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
   
   for(int i=0;i<1000;i++) {
      int a1 = i+5;
      int a2 = i*i;
      int res = func2(a1,a2);
      assert(res==a1+a2);
   }

   // -------------------------------- end tests

   #ifdef NDEBUG
   std::cout<< "WARNING: DEBUG DISABLED!\n";
   #else
   std::cout<< "TESTS COMPLETE.\n";
   #endif
}
