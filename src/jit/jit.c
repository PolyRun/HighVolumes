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
      jit_memoryPages_->pages = 100000;
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

jit_Table_8* jit_Table_8_prepend(jit_Table_8* old, uint8_t* bytes, uint8_t* src) {
   jit_Table_8* t = (jit_Table_8*)malloc(sizeof(jit_Table_8));
   if(old) {
      t->children = old->children+1;
   } else {
      t->children = 0;
   }
   t->next = old;
   t->src = src;
   for(int i=0;i<8;i++) {t->data[i] = bytes[i];}
   return t;
}

jit_Table_8* jit_immediate_8_via_data(const double val, const int xmm, jit_Table_8* t) {
   // c5 fb 10 05 xxxx  vmovsd 0x100(%rip),%xmm0
   // c5 fb 10 0d xxxx  vmovsd 0x100(%rip),%xmm1
   // c5 fb 10 15 xxxx  vmovsd 0x100(%rip),%xmm2
   // c5 fb 10 1d xxxx  vmovsd 0x100(%rip),%xmm3
   // c5 fb 10 25 xxxx  vmovsd 0x100(%rip),%xmm4
   // c5 fb 10 2d xxxx  vmovsd 0x100(%rip),%xmm5
   // c5 fb 10 35 xxxx  vmovsd 0x100(%rip),%xmm6
   // c5 fb 10 3d xxxx  vmovsd 0x100(%rip),%xmm7
   // c5 7b 10 05 xxxx  vmovsd 0x100(%rip),%xmm8
   // c5 7b 10 0d xxxx  vmovsd 0x100(%rip),%xmm9
   // c5 7b 10 15 xxxx  vmovsd 0x100(%rip),%xmm10
   // c5 7b 10 1d xxxx  vmovsd 0x100(%rip),%xmm11
   // c5 7b 10 25 xxxx  vmovsd 0x100(%rip),%xmm12
   // c5 7b 10 2d xxxx  vmovsd 0x100(%rip),%xmm13
   // c5 7b 10 35 xxxx  vmovsd 0x100(%rip),%xmm14
   // c5 7b 10 3d xxxx  vmovsd 0x100(%rip),%xmm15
   
   uint8_t b2 = 0xfb;
   if(xmm>7) {b2 = 0x7b;}
   uint8_t b4 = 0x05 + (xmm % 8)*8;
   {const uint8_t instr[] = {0xc5,b2,0x10,b4}; jit_push(instr,4);}
   // 32 bytes for the offset address
   {const uint8_t instr[] = {0xff,0xff,0xff,0xff}; jit_push(instr,4);}
   //{const uint8_t instr[] = {0,0,0,0}; jit_push(instr,4);}

   return jit_Table_8_prepend(t, (uint8_t*)&val, jit_head());
}

jit_Table_8* jit_broadcast_sd_via_data(const double val, const int ymm, jit_Table_8* t) {
// ASDF
   jit_vbroadcastsd_mem(jit_rip, 0xffffffff, ymm);
   return jit_Table_8_prepend(t, (uint8_t*)&val, jit_head());
}

void jit_table_8_consume(jit_Table_8* t) {
   if(t==NULL) {return;}

   jit_allign(8);
   int n = t->children+1;
   
   uint8_t* top = jit_head();
   double test = 1.0101010101;
   for(uint64_t i=0;i<n;i++) {// make space
      jit_push((const uint8_t*)&test,8);
   }
   
   int i = n; // write them in in reverse order, as list is prepend only
   while(t!=NULL) {
      i--;
      uint8_t* index = top+8*i;
      uint32_t offset = index - t->src;
      jit_write(t->src-4, (uint8_t*)&offset,4);
      jit_write(index, t->data, 8);

      jit_Table_8* next = t->next;
      free(t);
      t = next;
   }
   assert(i==0);
}

jit_Table_16* jit_Table_16_prepend(jit_Table_16* old, uint8_t* bytes, uint8_t* src) {
   jit_Table_16* t = (jit_Table_16*)malloc(sizeof(jit_Table_16));
   if(old) {
      t->children = old->children+1;
   } else {
      t->children = 0;
   }
   t->next = old;
   t->src = src;
   for(int i=0;i<16;i++) {t->data[i] = bytes[i];}
   return t;
}

jit_Table_16* jit_immediate_16_via_data(const double val0, const double val1, const int xmm, jit_Table_16* t) {
   // for 64 bytes:
   // c5 fd 28 05 00 01 00 	vmovapd 0x100(%rip),%ymm0        # 0x788
   // c5 fd 28 0d 00 01 00 	vmovapd 0x100(%rip),%ymm1        # 0x790
   // c5 fd 28 15 00 01 00 	vmovapd 0x100(%rip),%ymm2        # 0x798
   // c5 fd 28 1d 00 01 00 	vmovapd 0x100(%rip),%ymm3        # 0x7a0
   // c5 fd 28 25 00 01 00 	vmovapd 0x100(%rip),%ymm4        # 0x7a8
   // c5 fd 28 2d 00 01 00 	vmovapd 0x100(%rip),%ymm5        # 0x7b0
   // c5 fd 28 35 00 01 00 	vmovapd 0x100(%rip),%ymm6        # 0x7b8
   // c5 fd 28 3d 00 01 00 	vmovapd 0x100(%rip),%ymm7        # 0x7c0
   // c5 7d 28 05 00 01 00 	vmovapd 0x100(%rip),%ymm8        # 0x7c8
   // c5 7d 28 0d 00 01 00 	vmovapd 0x100(%rip),%ymm9        # 0x7d0
   // c5 7d 28 15 00 01 00 	vmovapd 0x100(%rip),%ymm10        # 0x7d8
   // c5 7d 28 1d 00 01 00 	vmovapd 0x100(%rip),%ymm11        # 0x7e0
   // c5 7d 28 25 00 01 00 	vmovapd 0x100(%rip),%ymm12        # 0x7e8
   // c5 7d 28 2d 00 01 00 	vmovapd 0x100(%rip),%ymm13        # 0x7f0
   // c5 7d 28 35 00 01 00 	vmovapd 0x100(%rip),%ymm14        # 0x7f8
   // c5 7d 28 3d 00 01 00 	vmovapd 0x100(%rip),%ymm15        # 0x800

   uint8_t b2 = 0xf9;
   if(xmm>7) {b2 = 0x79;}
   uint8_t b4 = 0x05 + (xmm % 8)*8;
   {const uint8_t instr[] = {0xc5,b2,0x28,b4}; jit_push(instr,4);}
   // 32 bytes for the offset address
   {const uint8_t instr[] = {0xff,0xff,0xff,0xff}; jit_push(instr,4);}
   //{const uint8_t instr[] = {0,0,0,0}; jit_push(instr,4);}
   double val[2] = {val0,val1};
   return jit_Table_16_prepend(t, (uint8_t*)&val, jit_head());
}

void jit_table_16_consume(jit_Table_16* t) {
   if(t==NULL) {return;}

   jit_allign(16);
   int n = t->children+1;
   
   uint8_t* top = jit_head();
   double test1 = 1.0101010101;
   double test2 = 2.0202020202;
   for(uint64_t i=0;i<n;i++) {// make space
      jit_push((const uint8_t*)&test1,8);
      jit_push((const uint8_t*)&test2,8);
   }
   
   int i = n; // write them in in reverse order, as list is prepend only
   while(t!=NULL) {
      i--;
      uint8_t* index = top+16*i;
      uint32_t offset = index - t->src;
      jit_write(t->src-4, (uint8_t*)&offset,4);
      jit_write(index, t->data, 16);

      jit_Table_16* next = t->next;
      free(t);
      t = next;
   }
   assert(i==0);
}

jit_Table_32* jit_Table_32_prepend(jit_Table_32* old, uint8_t* bytes, uint8_t* src) {
   jit_Table_32* t = (jit_Table_32*)malloc(sizeof(jit_Table_32));
   if(old) {
      t->children = old->children+1;
   } else {
      t->children = 0;
   }
   t->next = old;
   t->src = src;
   for(int i=0;i<32;i++) {t->data[i] = bytes[i];}
   return t;
}

jit_Table_32* jit_immediate_32_via_data(const double val0, const double val1, const double val2, const double val3, const int xmm, jit_Table_32* t) {
   // for 64 bytes:
   // c5 fd 28 05 00 01 00 	vmovapd 0x100(%rip),%ymm0        # 0x788
   // c5 fd 28 0d 00 01 00 	vmovapd 0x100(%rip),%ymm1        # 0x790
   // c5 fd 28 15 00 01 00 	vmovapd 0x100(%rip),%ymm2        # 0x798
   // c5 fd 28 1d 00 01 00 	vmovapd 0x100(%rip),%ymm3        # 0x7a0
   // c5 fd 28 25 00 01 00 	vmovapd 0x100(%rip),%ymm4        # 0x7a8
   // c5 fd 28 2d 00 01 00 	vmovapd 0x100(%rip),%ymm5        # 0x7b0
   // c5 fd 28 35 00 01 00 	vmovapd 0x100(%rip),%ymm6        # 0x7b8
   // c5 fd 28 3d 00 01 00 	vmovapd 0x100(%rip),%ymm7        # 0x7c0
   // c5 7d 28 05 00 01 00 	vmovapd 0x100(%rip),%ymm8        # 0x7c8
   // c5 7d 28 0d 00 01 00 	vmovapd 0x100(%rip),%ymm9        # 0x7d0
   // c5 7d 28 15 00 01 00 	vmovapd 0x100(%rip),%ymm10        # 0x7d8
   // c5 7d 28 1d 00 01 00 	vmovapd 0x100(%rip),%ymm11        # 0x7e0
   // c5 7d 28 25 00 01 00 	vmovapd 0x100(%rip),%ymm12        # 0x7e8
   // c5 7d 28 2d 00 01 00 	vmovapd 0x100(%rip),%ymm13        # 0x7f0
   // c5 7d 28 35 00 01 00 	vmovapd 0x100(%rip),%ymm14        # 0x7f8
   // c5 7d 28 3d 00 01 00 	vmovapd 0x100(%rip),%ymm15        # 0x800

   uint8_t b2 = 0xfd;
   if(xmm>7) {b2 = 0x7d;}
   uint8_t b4 = 0x05 + (xmm % 8)*8;
   {const uint8_t instr[] = {0xc5,b2,0x28,b4}; jit_push(instr,4);}
   // 32 bytes for the offset address
   {const uint8_t instr[] = {0xff,0xff,0xff,0xff}; jit_push(instr,4);}
   //{const uint8_t instr[] = {0,0,0,0}; jit_push(instr,4);}
   double val[4] = {val0,val1,val2,val3};
   return jit_Table_32_prepend(t, (uint8_t*)&val, jit_head());
}

void jit_table_32_consume(jit_Table_32* t) {
   if(t==NULL) {return;}

   jit_allign(32);
   int n = t->children+1;
   
   uint8_t* top = jit_head();
   double test1 = 1.0101010101;
   double test2 = 2.0202020202;
   double test3 = 3.0202020202;
   double test4 = 4.0202020202;
   for(uint64_t i=0;i<n;i++) {// make space
      jit_push((const uint8_t*)&test1,8);
      jit_push((const uint8_t*)&test2,8);
      jit_push((const uint8_t*)&test3,8);
      jit_push((const uint8_t*)&test4,8);
   }
   
   int i = n; // write them in in reverse order, as list is prepend only
   while(t!=NULL) {
      i--;
      uint8_t* index = top+32*i;
      uint32_t offset = index - t->src;
      jit_write(t->src-4, (uint8_t*)&offset,4);
      jit_write(index, t->data, 32);

      jit_Table_32* next = t->next;
      free(t);
      t = next;
   }
   assert(i==0);
}

void jit_permpd(uint8_t imm, int src, int dst) {
   uint8_t b2 = 0xe3;
   if(src < 8 && dst < 8) {
      b2 = 0xe3;
   } else if(src < 8 && dst >= 8) {
      b2 = 0x63;
   } else if(src >= 8 && dst < 8) {
      b2 = 0xc3;
   } else if(src >= 8 && dst >= 8) {
      b2 = 0x43;
   }
   uint8_t b5 = 0xc0 + (dst%8)*8 + (src%8)*1;
   
   // c4 e3 fd 01 c0 ff    	vpermpd $0xff,%ymm0,%ymm0
   { uint8_t instr[] = {0xc4,b2,0xfd,0x01,b5,imm}; jit_push(instr,6); }
}



void jit_permilpd(uint8_t imm, int src, int dst) {
   assert(imm<16);
   uint8_t b2 = 0xe3;
   if(src < 8 && dst < 8) {
      b2 = 0xe3;
   } else if(src < 8 && dst >= 8) {
      b2 = 0x63;
   } else if(src >= 8 && dst < 8) {
      b2 = 0xc3;
   } else if(src >= 8 && dst >= 8) {
      b2 = 0x43;
   }
   uint8_t b5 = 0xc0 + (dst%8)*8 + (src%8)*1;
   
   // c4 e3 79 05 c0 aa    	vpermilpd $0xaa,%xmm0,%xmm0
   { uint8_t instr[] = {0xc4,b2,0x7d,0x05,b5,imm}; jit_push(instr,6); }
}

void jit_permilpd_xmm(uint8_t imm, int src, int dst) {
   assert(imm<16);
   uint8_t b2 = 0xe3;
   if(src < 8 && dst < 8) {
      b2 = 0xe3;
   } else if(src < 8 && dst >= 8) {
      b2 = 0x63;
   } else if(src >= 8 && dst < 8) {
      b2 = 0xc3;
   } else if(src >= 8 && dst >= 8) {
      b2 = 0x43;
   }
   uint8_t b5 = 0xc0 + (dst%8)*8 + (src%8)*1;
   
   // c4 e3 79 05 c0 aa    	vpermilpd $0xaa,%xmm0,%xmm0
   { uint8_t instr[] = {0xc4,b2,0x79,0x05,b5,imm}; jit_push(instr,6); }
}


void jit_storeu_xmm(int src, jit_Register reg, uint32_t idx) {
   uint8_t b2 = 0xf9;
   if(src >= 8) {b2-=0x80;}
   uint8_t b4 = (src%8)*8;
   switch(reg) {
      case jit_rax: {b4+=0x80;break;}
      case jit_rcx: {b4+=0x81;break;}
      case jit_rdx: {b4+=0x82;break;}
      case jit_rbx: {b4+=0x83;break;}
      case jit_rsi: {b4+=0x86;break;}
      case jit_rdi: {b4+=0x87;break;}
      default: {assert(0 && "reg not handled!");}
   }
   // c5 f9 11 80 00 01 00 	vmovupd %xmm0,0x100(%rax)
   { uint8_t instr[] = {0xc5,b2,0x11,b4}; jit_push(instr,4); }
   jit_push((const uint8_t*)&idx,4);
}

void jit_storeu_ymm(int src, jit_Register reg, uint32_t idx) {
   uint8_t b2 = 0xfd;
   if(src >= 8) {b2-=0x80;}
   uint8_t b4 = (src%8)*8;
   switch(reg) {
      case jit_rax: {b4+=0x80;break;}
      case jit_rcx: {b4+=0x81;break;}
      case jit_rdx: {b4+=0x82;break;}
      case jit_rbx: {b4+=0x83;break;}
      case jit_rsi: {b4+=0x86;break;}
      case jit_rdi: {b4+=0x87;break;}
      default: {assert(0 && "reg not handled!");}
   }
   // c5 fd 11 80 00 01 00 	vmovupd %ymm0,0x100(%rax)
   { uint8_t instr[] = {0xc5,b2,0x11,b4}; jit_push(instr,4); }
   jit_push((const uint8_t*)&idx,4);
}

void jit_load_sd(jit_Register reg, uint32_t idx, int dst) {
   uint8_t b2 = 0xfb;
   if(dst >= 8) {b2-=0x80;}
   uint8_t b4 = (dst%8)*8;
   switch(reg) {
      case jit_rax: {b4+=0x80;break;}
      case jit_rcx: {b4+=0x81;break;}
      case jit_rdx: {b4+=0x82;break;}
      case jit_rbx: {b4+=0x83;break;}
      case jit_rsi: {b4+=0x86;break;}
      case jit_rdi: {b4+=0x87;break;}
      default: {assert(0 && "reg not handled!");}
   }
   // c5 fb 10 80 00 01 00 	vmovsd 0x100(%rax),%xmm0
   { uint8_t instr[] = {0xc5,b2,0x10,b4}; jit_push(instr,4); }
   jit_push((const uint8_t*)&idx,4);
}

void jit_loadu_xmm(jit_Register reg, uint32_t idx, int dst) {
   uint8_t b2 = 0xf9;
   if(dst >= 8) {b2-=0x80;}
   uint8_t b4 = (dst%8)*8;
   switch(reg) {
      case jit_rax: {b4+=0x80;break;}
      case jit_rcx: {b4+=0x81;break;}
      case jit_rdx: {b4+=0x82;break;}
      case jit_rbx: {b4+=0x83;break;}
      case jit_rsi: {b4+=0x86;break;}
      case jit_rdi: {b4+=0x87;break;}
      default: {assert(0 && "reg not handled!");}
   }
   // replace b3 with 0x28 to get aligned read
   // c5 f9 10 80 00 01 00 	vmovupd 0x100(%rax),%xmm0
   { uint8_t instr[] = {0xc5,b2,0x10,b4}; jit_push(instr,4); }
   jit_push((const uint8_t*)&idx,4);
}


void jit_loadu_ymm(jit_Register reg, uint32_t idx, int dst) {
   uint8_t b2 = 0xfd;
   if(dst >= 8) {b2-=0x80;}
   uint8_t b4 = (dst%8)*8;
   switch(reg) {
      case jit_rax: {b4+=0x80;break;}
      case jit_rcx: {b4+=0x81;break;}
      case jit_rdx: {b4+=0x82;break;}
      case jit_rbx: {b4+=0x83;break;}
      case jit_rsi: {b4+=0x86;break;}
      case jit_rdi: {b4+=0x87;break;}
      default: {assert(0 && "reg not handled!");}
   }
   // replace b3 with 0x28 to get aligned read
   // c5 fd 10 80 00 01 00 	vmovupd 0x100(%rax),%ymm0
   { uint8_t instr[] = {0xc5,b2,0x10,b4}; jit_push(instr,4); }
   jit_push((const uint8_t*)&idx,4);
}

void jit_vOPsd(uint8_t op, int src1, int src2, int dst) {
   // c5 fb 59 c0          	vmulsd %xmm0,%xmm0,%xmm0
   // c5 fb 5f c0          	vmaxsd %xmm0,%xmm0,%xmm0
   // c5 fb 5d c0          	vminsd %xmm0,%xmm0,%xmm0
   uint8_t b2 = 0xfb - src2*8;
   uint8_t b4 = 0xc0 + (src1%8)*1 + (dst%8)*8;
   if(dst >= 8) {b2-=0x80;}
   uint8_t bx = 0xc1;
   if(dst >= 8 && src1 >= 8) {bx-=0x80;}
   if(src1 < 8) {
      { uint8_t instr[] = {0xc5,b2,op,b4}; jit_push(instr,4); }
   } else {
      b2-=0x80;
      { uint8_t instr[] = {0xc4,bx,b2,op,b4}; jit_push(instr,5); }
   }
}

void jit_vmulsd(int src1, int src2, int dst) {
   // c5 f9 59 c0          	vmulpd %xmm0,%xmm0,%xmm0
   jit_vOPsd(0x59,src1,src2,dst);
}

void jit_vmaxsd(int src1, int src2, int dst) {
   // c5 f9 5f c0          	vmaxpd %xmm0,%xmm0,%xmm0
   jit_vOPsd(0x5f,src1,src2,dst);
}

void jit_vminsd(int src1, int src2, int dst) {
   jit_vOPsd(0x5d,src1,src2,dst);
}


void jit_vOPsd_mem(uint8_t op, jit_Register reg, uint32_t idx, int src2, int dst) {
   // c5 fb 59 80 00 01 00 	vmulsd 0x100(%rax),%xmm0,%xmm0
   // c5 f9 59 80 00 01 00 	vmulpd 0x100(%rax),%xmm0,%xmm0
   uint8_t b2 = 0xfb - src2*8;
   uint8_t b4 = (dst%8)*8;
   if(dst >= 8) {b2-=0x80;}
   
   switch(reg) {
      case jit_rax: {b4+=0x80;break;}
      case jit_rcx: {b4+=0x81;break;}
      case jit_rdx: {b4+=0x82;break;}
      case jit_rbx: {b4+=0x83;break;}
      case jit_rsi: {b4+=0x86;break;}
      case jit_rdi: {b4+=0x87;break;}
      default: {assert(0 && "reg not handled!");}
   }
   
   { uint8_t instr[] = {0xc5,b2,op,b4}; jit_push(instr,4); }
   jit_push((const uint8_t*)&idx,4);
}

void jit_vmulsd_mem(jit_Register reg, uint32_t idx, int src2, int dst) {
   // c5 f9 59 80 00 01 00 	vmulpd 0x100(%rax),%xmm0,%xmm0
   jit_vOPsd_mem(0x59,reg,idx,src2,dst);
}


void jit_vOPpd_xmm(uint8_t op, int src1, int src2, int dst) {
   // c5 f9 59 c0          	vmulpd %xmm0,%xmm0,%xmm0
   uint8_t b2 = 0xf9 - src2*8;
   uint8_t b4 = 0xc0 + (src1%8)*1 + (dst%8)*8;
   if(dst >= 8) {b2-=0x80;}
   uint8_t bx = 0xc1;
   if(dst >= 8 && src1 >= 8) {bx-=0x80;}
   if(src1 < 8) {
      { uint8_t instr[] = {0xc5,b2,op,b4}; jit_push(instr,4); }
   } else {
      b2-=0x80;
      { uint8_t instr[] = {0xc4,bx,b2,op,b4}; jit_push(instr,5); }
   }
}

void jit_vmulpd_xmm(int src1, int src2, int dst) {
   // c5 f9 59 c0          	vmulpd %xmm0,%xmm0,%xmm0
   jit_vOPpd_xmm(0x59,src1,src2,dst);
}

void jit_vmaxpd_xmm(int src1, int src2, int dst) {
   // c5 f9 5f c0          	vmaxpd %xmm0,%xmm0,%xmm0
   jit_vOPpd_xmm(0x5f,src1,src2,dst);
}

void jit_vminpd_xmm(int src1, int src2, int dst) {
   jit_vOPpd_xmm(0x5d,src1,src2,dst);
}

void jit_vOPpd_mem_xmm(uint8_t op, jit_Register reg, uint32_t idx, int src2, int dst) {
   // c5 f9 59 c0          	vmulpd %xmm0,%xmm0,%xmm0
   // c5 f9 59 80 00 01 00 	vmulpd 0x100(%rax),%xmm0,%xmm0
   uint8_t b2 = 0xf9 - src2*8;
   uint8_t b4 = (dst%8)*8;
   if(dst >= 8) {b2-=0x80;}
   
   switch(reg) {
      case jit_rax: {b4+=0x80;break;}
      case jit_rcx: {b4+=0x81;break;}
      case jit_rdx: {b4+=0x82;break;}
      case jit_rbx: {b4+=0x83;break;}
      case jit_rsi: {b4+=0x86;break;}
      case jit_rdi: {b4+=0x87;break;}
      default: {assert(0 && "reg not handled!");}
   }
   
   { uint8_t instr[] = {0xc5,b2,op,b4}; jit_push(instr,4); }
   jit_push((const uint8_t*)&idx,4);
}

void jit_fma_mem(jit_Register reg, uint32_t idx, int a, int b, uint8_t b3, uint8_t b4) {
   // c4 e2 f9 a9 80 00 01 	vfmadd213sd 0x100(%rax),%xmm0,%xmm0
   // c4 e2 f9 a8 80 00 01 	vfmadd213pd 0x100(%rax),%xmm0,%xmm0
   // c4 e2 fd a8 80 00 01 	vfmadd213pd 0x100(%rax),%ymm0,%ymm0
   // c4 e2 fd 98 80 00 01 	vfmadd132pd 0x100(%rax),%ymm0,%ymm0
   // c4 e2 fd b8 80 00 01 	vfmadd231pd 0x100(%rax),%ymm0,%ymm0
   uint8_t b5 = 8*(b%8);
   uint8_t b2 = 0xe2;
   if(b>=8) {b2-=0x80;}
   b3-= a*8;
   switch(reg) {
      case jit_rax: {b5+=0x80;break;}
      case jit_rcx: {b5+=0x81;break;}
      case jit_rdx: {b5+=0x82;break;}
      case jit_rbx: {b5+=0x83;break;}
      case jit_rsi: {b5+=0x86;break;}
      case jit_rdi: {b5+=0x87;break;}
      default: {assert(0 && "reg not handled!");}
   }
   { uint8_t instr[] = {0xc4,b2,b3,b4,b5}; jit_push(instr,5); }
   jit_push((const uint8_t*)&idx,4);
}

void jit_vfmad213sd_mem(jit_Register reg, uint32_t idx, int src, int dst) {
   // c4 e2 f9 a9 80 00 01 	vfmadd213sd 0x100(%rax),%xmm0,%xmm0
   jit_fma_mem(reg,idx,src,dst,0xf9,0xa9);
}

void jit_vfmad213pd_mem_xmm(jit_Register reg, uint32_t idx, int src, int dst) {
   // c4 e2 f9 a8 80 00 01 	vfmadd213pd 0x100(%rax),%xmm0,%xmm0
   jit_fma_mem(reg,idx,src,dst,0xf9,0xa8);
}

void jit_vfmad213pd_mem_ymm(jit_Register reg, uint32_t idx, int src, int dst) {
   // c4 e2 fd a8 80 00 01 	vfmadd213pd 0x100(%rax),%ymm0,%ymm0
   jit_fma_mem(reg,idx,src,dst,0xfd,0xa8);
}

void jit_vfnmad231pd_mem_ymm(jit_Register reg, uint32_t idx, int src, int dst) {
   // c4 e2 f9 bc 80 00 01 	vfnmadd231pd 0x100(%rax),%xmm0,%xmm0
   // 	c4 e2 fd bc 80 00 01 	vfnmadd231pd 0x100(%rax),%ymm0,%ymm0
   jit_fma_mem(reg,idx,src,dst,0xfd,0xbc);
}

void jit_vmulpd_mem_xmm(jit_Register reg, uint32_t idx, int src2, int dst) {
   // c5 f9 59 80 00 01 00 	vmulpd 0x100(%rax),%xmm0,%xmm0
   jit_vOPpd_mem_xmm(0x59,reg,idx,src2,dst);
}

void jit_vOPpd_ymm(uint8_t op, int src1, int src2, int dst) {
   // c5 f9 59 c0          	vmulpd %xmm0,%xmm0,%xmm0
   uint8_t b2 = 0xfd - src2*8;
   uint8_t b4 = 0xc0 + (src1%8)*1 + (dst%8)*8;
   if(dst >= 8) {b2-=0x80;}
   uint8_t bx = 0xc1;
   if(dst >= 8 && src1 >= 8) {bx-=0x80;}
   if(src1 < 8) {
      { uint8_t instr[] = {0xc5,b2,op,b4}; jit_push(instr,4); }
   } else {
      b2-=0x80;
      { uint8_t instr[] = {0xc4,bx,b2,op,b4}; jit_push(instr,5); }
   }
}

void jit_vmulpd_ymm(int src1, int src2, int dst) {
   // c5 f9 59 c0          	vmulpd %xmm0,%xmm0,%xmm0
   jit_vOPpd_ymm(0x59,src1,src2,dst);
}

void jit_vmaxpd_ymm(int src1, int src2, int dst) {
   // c5 f9 5f c0          	vmaxpd %xmm0,%xmm0,%xmm0
   jit_vOPpd_ymm(0x5f,src1,src2,dst);
}

void jit_vminpd_ymm(int src1, int src2, int dst) {
   jit_vOPpd_ymm(0x5d,src1,src2,dst);
}

void jit_vOPpd_mem_ymm(uint8_t op, jit_Register reg, uint32_t idx, int src2, int dst) {
   // c5 f9 59 c0          	vmulpd %xmm0,%xmm0,%xmm0
   // c5 f9 59 80 00 01 00 	vmulpd 0x100(%rax),%xmm0,%xmm0
   // c5 fd 59 80 00 01 00 	vmulpd 0x100(%rax),%ymm0,%ymm0
   uint8_t b2 = 0xfd - src2*8;
   uint8_t b4 = (dst%8)*8;
   if(dst >= 8) {b2-=0x80;}
   
   switch(reg) {
      case jit_rax: {b4+=0x80;break;}
      case jit_rcx: {b4+=0x81;break;}
      case jit_rdx: {b4+=0x82;break;}
      case jit_rbx: {b4+=0x83;break;}
      case jit_rsi: {b4+=0x86;break;}
      case jit_rdi: {b4+=0x87;break;}
      default: {assert(0 && "reg not handled!");}
   }
   
   { uint8_t instr[] = {0xc5,b2,op,b4}; jit_push(instr,4); }
   jit_push((const uint8_t*)&idx,4);
}



void jit_vmulpd_mem_ymm(jit_Register reg, uint32_t idx, int src2, int dst) {
   // c5 f9 59 80 00 01 00 	vmulpd 0x100(%rax),%xmm0,%xmm0
   // c5 fd 59 80 00 01 00 	vmulpd 0x100(%rax),%ymm0,%ymm0
   jit_vOPpd_mem_ymm(0x59,reg,idx,src2,dst);
}

void jit_vbroadcastsd_ymm(int src, int dst) {
   // c4 e2 7d 19 c0       	vbroadcastsd %xmm0,%ymm0
   uint8_t b2 = 0xe2;
   uint8_t b5 = 0xc0 + 8*(dst%8) + 1*(src%8);
   if(dst>=8) {b2-=0x80;}
   if(src>=8) {b2-=0x20;}
   { uint8_t instr[] = {0xc4,b2,0x7d,0x19,b5}; jit_push(instr,5); }
}

void jit_vbroadcastsd_mem(jit_Register reg, uint32_t idx, int dst) {
   uint8_t b2 = 0xe2;
   uint8_t b5 = (dst%8)*8;
   if(dst >= 8) {b2-=0x80;}
   switch(reg) {
      case jit_rax: {b5+=0x80;break;}
      case jit_rcx: {b5+=0x81;break;}
      case jit_rdx: {b5+=0x82;break;}
      case jit_rbx: {b5+=0x83;break;}
      case jit_rsi: {b5+=0x86;break;}
      case jit_rdi: {b5+=0x87;break;}
      case jit_rip: {b5+=0x05;break;}
      default: {assert(0 && "reg not handled!");}
   }
   
   // c4 e2 7d 19 80 00 01 	vbroadcastsd 0x100(%rax),%ymm0
   { uint8_t instr[] = {0xc4,b2,0x7d,0x19,b5}; jit_push(instr,5); }
   jit_push((const uint8_t*)&idx,4);
}

void jit_emit_vzeroupper() {
   // c5 f8 77                vzeroupper
   { uint8_t instr[] = {0xc5,0xf8,0x77}; jit_push(instr,3); }
}


void jit_emit_return() {
   { uint8_t instr[] = {0xf3,0xc3}; jit_push(instr,2); }
}

