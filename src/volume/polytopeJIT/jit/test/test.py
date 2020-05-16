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
 



