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
#include <assert.h>

#define JIT_RDI 


typedef struct jit_MemoryPages {
   uint8_t *mem;      // pointer to memory region
   size_t page_size;  // size of mem pages
   size_t pages;      // number of pages we have at mem
   size_t head;       // offset to mem, position of first non-used
} jit_MemoryPages;

jit_MemoryPages* jit_memoryPages();

// calculate how many pages one needs to fit bytes
int jit_estimate_npages(const size_t bytes_requested);

// push character to mem
void jit_pushByte(const uint8_t c);
// push array of n chars to mem
void jit_push(const uint8_t* c, const size_t n);
// write n bytes from c to mem
void jit_write(uint8_t* mem, const uint8_t* c, const size_t n);

// add zero, until alligned
void jit_allign(const size_t a);

// print mem to stdout
void jit_print();

uint8_t* jit_head();

#endif // HEADER_JIT_H
