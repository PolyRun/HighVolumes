#!/usr/bin/env python3


for i in range(16):
    for j in range(16):
        print("vpermilpd $0b10101010,%ymm{},%ymm{}".format(i,j))

for i in range(16):
    print("vmovapd 0x100(%rip),%xmm{}".format(i));

for i in range(16):
    print("vmovapd 0x100(%rip),%ymm{}".format(i));
