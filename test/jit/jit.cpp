#include <iostream>

extern "C" { // must be included C stlye
#include "../../src/jit/jit.h"
}

int main() {
   std::cout << "start test.\n";

   // -------------------------------- start tests
   printf("memory: %d\n",jit_estimate_memory(1));
   printf("memory: %d\n",jit_estimate_memory(4096 + 1));
   // -------------------------------- end tests

   #ifdef NDEBUG
   std::cout<< "WARNING: DEBUG DISABLED!\n";
   #else
   std::cout<< "TESTS COMPLETE.\n";
   #endif
}
