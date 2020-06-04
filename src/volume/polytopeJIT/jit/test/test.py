#!/usr/bin/env python3


for i in range(16):
    for j in range(16):
        print("vpermilpd $0b10101010,%ymm{},%ymm{}".format(i,j))

for i in range(16):
    print("vmovapd 0x100(%rip),%xmm{}".format(i));

for i in range(16):
    print("vmovapd 0x100(%rip),%ymm{}".format(i));

for i in range(16):
    print("vmovupd 0x100(%rax),%xmm{}".format(i));
    print("vmovupd 0x100(%rbx),%xmm{}".format(i));
    print("vmovupd 0x100(%rcx),%xmm{}".format(i));
    print("vmovupd 0x100(%rdx),%xmm{}".format(i));
    print("vmovupd 0x100(%rdi),%xmm{}".format(i));
    print("vmovupd 0x100(%rsi),%xmm{}".format(i));

print("########################");
for i in range(16):
    for j in range(16):
        print("vmulpd %xmm{},%xmm{},%xmm8".format(i,j))
 
print("vmulpd %xmm2,%xmm2,%xmm2")

print("vmulpd %xmm2,%xmm3,%xmm7")
print("vmulpd %xmm2,%xmm3,%xmm8")

print("vmulpd %xmm8,%xmm2,%xmm7")
print("vmulpd %xmm8,%xmm2,%xmm8")

for i in range(16):
    for j in range(16):
        print("vmaxpd %xmm{},%xmm{},%xmm0".format(i,j))

print("vmulpd %xmm0,%xmm0,%xmm0".format(i,j))
print("vmaxpd %xmm0,%xmm0,%xmm0".format(i,j))
print("vminpd %xmm0,%xmm0,%xmm0".format(i,j))
 
for i in range(16):
    print("vmovupd 0x100(%rax),%ymm{}".format(i));

print("vmulpd %xmm0,%xmm0,%xmm0".format(i,j))
print("vmaxpd %xmm0,%xmm0,%xmm0".format(i,j))
print("vminpd %xmm0,%xmm0,%xmm0".format(i,j))
 
print("vmulpd %ymm0,%ymm0,%ymm0".format(i,j))
print("vmaxpd %ymm0,%ymm0,%ymm0".format(i,j))
print("vminpd %ymm0,%ymm0,%ymm0".format(i,j))
 

print("########################");
for i in range(16):
    for j in range(16):
        print("vmulpd 0x100(%rax),%xmm{},%xmm{}".format(i,j))
        print("vmulpd 0x100(%rbx),%xmm{},%xmm{}".format(i,j))
        print("vmulpd 0x100(%rcx),%xmm{},%xmm{}".format(i,j))
        print("vmulpd 0x100(%rdx),%xmm{},%xmm{}".format(i,j))
        print("vmulpd 0x100(%rdi),%xmm{},%xmm{}".format(i,j))
        print("vmulpd 0x100(%rsi),%xmm{},%xmm{}".format(i,j))


for i in range(16):
    for j in range(16):
        print("vpermilpd $0b10101010,%xmm{},%xmm{}".format(i,j))
        print("vpermilpd $0b10101010,%ymm{},%ymm{}".format(i,j))



print("vmulpd %xmm0,%xmm0,%xmm0".format(i,j))
print("vmaxpd %xmm0,%xmm0,%xmm0".format(i,j))
print("vminpd %xmm0,%xmm0,%xmm0".format(i,j))
 
print("vmulsd %xmm0,%xmm0,%xmm0".format(i,j))
print("vmaxsd %xmm0,%xmm0,%xmm0".format(i,j))
print("vminsd %xmm0,%xmm0,%xmm0".format(i,j))
 

print("vmulsd 0x100(%rax),%xmm{},%xmm{}".format(0,0))
print("vmulpd 0x100(%rax),%xmm{},%xmm{}".format(0,0))

print("vextractf128 $0,%ymm0,%xmm0")
print("vperm2f128 $0b00010001,%ymm0,%ymm0,%ymm0")

#for i in range(16):
#    for j in range(16):
#       print("vperm2f128 $0b00010001,%ymm0,%ymm0,%ymm0")

print("vpermpd $0b11111111,%ymm0,%ymm0")

print("vmulpd 0x100(%rax),%ymm0,%ymm0")

for i in range(16):
   print("vmovupd %ymm{},0x100(%rax)".format(i))
   print("vmovupd %xmm{},0x100(%rax)".format(i))
   print("vmovupd %xmm{},0x100(%rcx)".format(i))



print("########################");
for i in range(16):
    for j in range(16):
        print("vfmadd213sd 0x100(%rax),%xmm{},%xmm{}".format(i,j))
        #print("vfmadd213pd 0x100(%rax),%xmm{},%xmm{}".format(i,j))
        #print("vfmadd213pd 0x100(%rax),%ymm{},%ymm{}".format(i,j))
        #print("vfmadd132pd 0x100(%rax),%ymm{},%ymm{}".format(i,j))
        #print("vfmadd231pd 0x100(%rax),%ymm{},%ymm{}".format(i,j))
        
        #print("vmulpd 0x100(%rax),%xmm{},%xmm{}".format(i,j))
        #print("vmulpd 0x100(%rbx),%xmm{},%xmm{}".format(i,j))
        #print("vmulpd 0x100(%rcx),%xmm{},%xmm{}".format(i,j))
        #print("vmulpd 0x100(%rdx),%xmm{},%xmm{}".format(i,j))
        #print("vmulpd 0x100(%rdi),%xmm{},%xmm{}".format(i,j))
        #print("vmulpd 0x100(%rsi),%xmm{},%xmm{}".format(i,j))


print("########################");
for i in range(16):
    for j in range(16):
        print("vbroadcastsd %xmm{},%ymm{}".format(i,j))


for i in range(16):
    print("vmovsd 0x100(%rax),%xmm{}".format(i));
    print("vmovsd 0x100(%rbx),%xmm{}".format(i));
    print("vmovsd 0x100(%rcx),%xmm{}".format(i));
    print("vmovsd 0x100(%rdx),%xmm{}".format(i));
    print("vmovsd 0x100(%rdi),%xmm{}".format(i));
    print("vmovsd 0x100(%rsi),%xmm{}".format(i));

for i in range(16):
    print("vbroadcastsd 0x100(%rax), %ymm{}".format(i));
    print("vbroadcastsd 0x100(%rbx), %ymm{}".format(i));
    print("vbroadcastsd 0x100(%rip), %ymm{}".format(i));


print("########################");
for i in range(16):
    for j in range(16):
        print("vfnmadd231pd 0x100(%rax),%ymm{},%ymm{}".format(i,j))



