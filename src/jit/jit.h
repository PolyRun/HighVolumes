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

// resets the head to mem.
// invalidates all functions
void jit_clear();

// ---------------------------------- op generators
void jit_immediate_via_rax(const double val, const int reg);

// table: 1 double
typedef struct jit_Table_8 {
   struct jit_Table_8* next; // linked list
   size_t children; // number of children in list
   uint8_t data[8]; // payload
   uint8_t* src; // end of instruction
   // will write relative offset at src-4;
} jit_Table_8;

jit_Table_8* jit_Table_8_prepend(jit_Table_8* old, uint8_t* bytes, uint8_t* src);

// sets op down
// prepends entry to table
jit_Table_8* jit_immediate_via_data(const double val, const int xmm, jit_Table_8* t);

// consume and free table
// write table in memory, go set up references to it
void jit_table_consume(jit_Table_8* t);

void jit_emit_return();

#endif // HEADER_JIT_H
