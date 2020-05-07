#include <iostream>
#include <string>

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
   jit_push(0x48); // mov $1, %rax # set for syscall write
   jit_push(0xc7);
   jit_push(0xc0);
   jit_push(0x01);
   jit_push(0x00);
   jit_push(0x00);
   jit_push(0x00);
   #elif __APPLE__
   jit_push(0x48);
   jit_push(0xc7);
   jit_push(0xc0);
   jit_push(0x04);
   jit_push(0x00);
   jit_push(0x00);
   jit_push(0x02);
   #endif
  
   jit_push(0x48); // mov $1,rdi   # to stdout
   jit_push(0xc7);
   jit_push(0xc7);
   jit_push(0x01);
   jit_push(0x00);
   jit_push(0x00);
   jit_push(0x00);

   jit_push(0x48); // offset to end of code / start of string
   jit_push(0x8d);
   jit_push(0x35);
   jit_push(0x0a); // this is the offset, starts at end of this instr
   jit_push(0x00);
   jit_push(0x00);
   jit_push(0x00);
   
   jit_push(0x48); // mov msg_size, %rdx  # store message size
   jit_push(0xc7);
   jit_push(0xc2);
   jit_push((msg_size & 0xFF) >> 0);// set this
   jit_push((msg_size & 0xFF00) >> 8);// set this
   jit_push((msg_size & 0xFF0000) >> 16);// set this
   jit_push((msg_size & 0xFF000000) >> 24);// set this
   
   jit_push(0x0f); // syscall
   jit_push(0x05);
 
   jit_push(0xc3);  // ret
   
   for(auto c : msg) { // add string data at end
      jit_push(c);
   }

   jit_print();

   func(); // call the function!
   // -------------------------------- end tests

   #ifdef NDEBUG
   std::cout<< "WARNING: DEBUG DISABLED!\n";
   #else
   std::cout<< "TESTS COMPLETE.\n";
   #endif
}
