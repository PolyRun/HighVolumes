#include "inside.h"

/*
# compile c to asm:
# gcc -S -fverbose-asm -O2 foo.c

# compile asm to machinecode:
# gcc -c hello.s

# dump machinecode to hex: 
# objdump -D hello.o

find test scripts in test/

*/

void PolytopeJIT_generate_inside_ref(const Polytope *p, PolytopeJIT *o) {
   jit_print();
   
   // info about calling convention:
   // input pointer -> %rdi
   // output boolean -> %eax
   
   // required computation:
   // for each constraint i:
   //    check ai*x <= b

   o->inside = (pjit_inside_f_t)jit_head();
   
   // return false;
   // ------------------
   // xorl %eax, %eax
   // ret
   // ------------------
   // 31 c0
   // c3
   { uint8_t instr[] = {0x31,0xc0,0xc3}; jit_push(instr,3); }
   jit_print();
   // TODO: implement correctly!
}

