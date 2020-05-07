#include "jit.h"

jit_MemoryPages* jit_memoryPages = NULL;

int jit_estimate_memory(const size_t bytes_requested) {
   size_t psize = sysconf(_SC_PAGE_SIZE);
   size_t npages = (bytes_requested + psize - 1)/psize;
   return npages * psize;
}

