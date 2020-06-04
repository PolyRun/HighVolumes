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

// ------------------- registers
typedef enum jit_Register {
   jit_rdi,// arg 1
   jit_rsi,// arg 2
   jit_rdx,// arg 3
   jit_rcx,// arg 4
   jit_rax,
   jit_rbx,
   jit_rip,
} jit_Register;

// -------------------------------------- table 8
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
jit_Table_8* jit_immediate_8_via_data(const double val, const int xmm, jit_Table_8* t);
jit_Table_8* jit_broadcast_sd_via_data(const double val, const int ymm, jit_Table_8* t);

// consume and free table
// write table in memory, go set up references to it
void jit_table_8_consume(jit_Table_8* t);

// -------------------------------------- table 16
// table: 2 double
typedef struct jit_Table_16 {
   struct jit_Table_16* next; // linked list
   size_t children;  // number of children in list
   uint8_t data[16]; // payload
   uint8_t* src;     // end of instruction
   // will write relative offset at src-4;
} jit_Table_16;

jit_Table_16* jit_Table_16_prepend(jit_Table_16* old, uint8_t* bytes, uint8_t* src);

// sets op down
// prepends entry to table
jit_Table_16* jit_immediate_16_via_data(const double val0, const double val1, const int xmm, jit_Table_16* t);

// consume and free table
// write table in memory, go set up references to it
void jit_table_16_consume(jit_Table_16* t);

// -------------------------------------- table 32
// table: 4 double
typedef struct jit_Table_32 {
   struct jit_Table_32* next; // linked list
   size_t children;  // number of children in list
   uint8_t data[32]; // payload
   uint8_t* src;     // end of instruction
   // will write relative offset at src-4;
} jit_Table_32;

jit_Table_32* jit_Table_32_prepend(jit_Table_32* old, uint8_t* bytes, uint8_t* src);

// sets op down
// prepends entry to table
jit_Table_32* jit_immediate_32_via_data(const double val0, const double val1, const double val2, const double val3, const int xmm, jit_Table_32* t);

// consume and free table
// write table in memory, go set up references to it
void jit_table_32_consume(jit_Table_32* t);



//------------------------------------- flops:
// 3 lat, _mm256_permute4x64_pd
void jit_permpd(uint8_t imm, int src, int dst);
// 1 lat, _mm256_permute_pd
void jit_permilpd(uint8_t imm, int src, int dst);
void jit_permilpd_xmm(uint8_t imm, int src, int dst);

void jit_load_sd(jit_Register reg, uint32_t idx, int dst);
void jit_vmulsd(int src1, int src2, int dst);
void jit_vmulsd_mem(jit_Register reg, uint32_t idx, int src2, int dst);
void jit_vmaxsd(int src1, int src2, int dst);
void jit_vminsd(int src1, int src2, int dst);

void jit_loadu_xmm(jit_Register reg, uint32_t idx, int dst);
void jit_storeu_xmm(int src, jit_Register reg, uint32_t idx);
void jit_vmulpd_xmm(int src1, int src2, int dst);
void jit_vmulpd_mem_xmm(jit_Register reg, uint32_t idx, int src2, int dst);
void jit_vmaxpd_xmm(int src1, int src2, int dst);
void jit_vminpd_xmm(int src1, int src2, int dst);

void jit_loadu_ymm(jit_Register reg, uint32_t idx, int dst);
void jit_storeu_ymm(int src, jit_Register reg, uint32_t idx);
void jit_vmulpd_ymm(int src1, int src2, int dst);
void jit_vmulpd_mem_ymm(jit_Register reg, uint32_t idx, int src2, int dst);
void jit_vmaxpd_ymm(int src1, int src2, int dst);
void jit_vminpd_ymm(int src1, int src2, int dst);

// dst = idx(%reg) + dst*src;
void jit_vfmad213sd_mem(jit_Register reg, uint32_t idx, int src, int dst);
void jit_vfmad213pd_mem_xmm(jit_Register reg, uint32_t idx, int src, int dst);
void jit_vfmad213pd_mem_ymm(jit_Register reg, uint32_t idx, int src, int dst);

// dst = dst - idx(%reg)*src
void jit_vfnmad231pd_mem_ymm(jit_Register reg, uint32_t idx, int src, int dst);

void jit_vbroadcastsd_ymm(int src, int dst);
void jit_vbroadcastsd_mem(jit_Register reg, uint32_t idx, int dst);

void jit_emit_vzeroupper();

void jit_emit_return();

#endif // HEADER_JIT_H
