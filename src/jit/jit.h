// Here we are building a facility to push machine code to executable memory.
// inspiration taken from:
// https://solarianprogrammer.com/2018/01/12/writing-minimal-x86-64-jit-compiler-cpp-part-2/

#ifndef HEADER_JIT_H
#define HEADER_JIT_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/mman.h>
#include <stdint.h>

typedef struct jit_MemoryPages {
   uint8_t *mem;      // pointer to memory region
   size_t page_size;  // size of mem pages
   size_t pages;      // number of pages we have at mem
   size_t head;       // offset to mem, position of first non-used
} jit_MemoryPages;

extern jit_MemoryPages* jit_memoryPages;

int jit_estimate_memory(const size_t bytes_requested);

#endif // HEADER_JIT_H
