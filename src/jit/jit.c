#include "jit.h"

jit_MemoryPages* jit_memoryPages_ = NULL;

uint8_t* jit_mmap(const size_t npages) {
   size_t psize = sysconf(_SC_PAGE_SIZE);
   uint8_t *mem = (uint8_t*) mmap(NULL, npages*psize, PROT_READ | PROT_WRITE | PROT_EXEC, MAP_PRIVATE | MAP_ANONYMOUS ,-1, 0);
   return mem;
}


jit_MemoryPages* jit_memoryPages() {
   if(!jit_memoryPages_) {
      jit_memoryPages_ = (jit_MemoryPages*)malloc(sizeof(jit_MemoryPages));
      
      size_t psize = sysconf(_SC_PAGE_SIZE);
      jit_memoryPages_->page_size = psize;
      jit_memoryPages_->pages = 500;
      jit_memoryPages_->head = 0;
      jit_memoryPages_->mem = jit_mmap(jit_memoryPages_->pages);;
   }

   assert(jit_memoryPages_->head < jit_memoryPages_->page_size * jit_memoryPages_->pages && "out of mmap memory!");
   
   return jit_memoryPages_;
}

int jit_estimate_npages(const size_t bytes_requested) {
   size_t psize = sysconf(_SC_PAGE_SIZE);
   return (bytes_requested + psize - 1)/psize;
}

void jit_pushByte(const uint8_t c) {
   jit_MemoryPages* p = jit_memoryPages();
   p->mem[p->head] = c;
   p->head++;
}

void jit_push(const uint8_t* c, const size_t n) {
   jit_MemoryPages* p = jit_memoryPages();
   
   for(int i=0;i<n;i++) {
      assert(jit_memoryPages_->head < jit_memoryPages_->page_size * jit_memoryPages_->pages && "out of mmap memory!");
      p->mem[p->head] = c[i];
      p->head++;
   }
}

void jit_write(uint8_t* mem, const uint8_t* c, const size_t n) {
   for(int i=0;i<n;i++) {
      mem[i] = c[i];
   }
}

void jit_allign(const size_t a) {
   jit_MemoryPages* p = jit_memoryPages();
   while(p->head % a != 0) {jit_pushByte(0);}
}

void jit_print() {
   jit_MemoryPages* p = jit_memoryPages();
   printf("Memory content: %p %ld\n", p->mem, p->head);
   for(int i=0;i<p->head;i++) {
      uint8_t c = p->mem[i];
      printf("%x ",c);
      if(i%16==0 && i>0) {
         printf("\n");
      }
   }
   printf("\n");
}

uint8_t* jit_head() {
   jit_MemoryPages* p = jit_memoryPages();
   return p->mem + p->head;
}


void jit_clear() {
   jit_MemoryPages* p = jit_memoryPages();
   p->head = 0;
}

void jit_immediate_via_rax(const double val, const int xmm) {
   //movabs $0xff00ff00ff00ff00,%rax
   {const uint8_t instr[] = {0x48,0xb8}; jit_push(instr,2); }
   jit_push((const uint8_t*)&val,8);
   // c4 e1 f9 6e c0       	vmovq  %rax,%xmm0
   // c4 e1 f9 6e c8       	vmovq  %rax,%xmm1
   // c4 e1 f9 6e d0       	vmovq  %rax,%xmm2
   // c4 e1 f9 6e d8       	vmovq  %rax,%xmm3
   // c4 e1 f9 6e e0       	vmovq  %rax,%xmm4
   // c4 e1 f9 6e e8       	vmovq  %rax,%xmm5
   // c4 e1 f9 6e f0       	vmovq  %rax,%xmm6
   // c4 e1 f9 6e f8       	vmovq  %rax,%xmm7
   // c4 61 f9 6e c0       	vmovq  %rax,%xmm8
   // c4 61 f9 6e c8       	vmovq  %rax,%xmm9
   // c4 61 f9 6e d0       	vmovq  %rax,%xmm10
   // c4 61 f9 6e d8       	vmovq  %rax,%xmm11
   // c4 61 f9 6e e0       	vmovq  %rax,%xmm12
   // c4 61 f9 6e e8       	vmovq  %rax,%xmm13
   // c4 61 f9 6e f0       	vmovq  %rax,%xmm14
   // c4 61 f9 6e f8       	vmovq  %rax,%xmm15

   assert(xmm < 16 && xmm >= 0);
   uint8_t b2 = 0xe1;
   if(xmm > 7) {b2 = 0x61;}
   uint8_t b5 = 0xc0 + (xmm % 8)*8;
   {const uint8_t instr[] = {0xc4,b2,0xf9,0x6e,b5}; jit_push(instr,5);}
}


