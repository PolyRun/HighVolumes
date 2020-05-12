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
      jit_memoryPages_->pages = 100;
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


