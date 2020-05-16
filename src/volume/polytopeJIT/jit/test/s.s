.global _inside

.text
_inside:
xorl %eax, %eax
vmovq  %rsi,%xmm0
vmovq  %rsi,%xmm1
vmovq  %rsi,%xmm2
ja L_end
jmp L_end
mulsd   8(%rdi), %xmm0
mulsd   16(%rdi), %xmm0
mulsd   512(%rdi), %xmm0
mulsd   24(%rdi), %xmm1
mulsd   512(%rdi), %xmm1
mulsd   512(%rdi), %xmm1
mulsd   512(%rdi), %xmm1
mulsd   512(%rdi), %xmm1
mulsd   512(%rdi), %xmm1
mulsd   512(%rdi), %xmm1
mulsd   512(%rdi), %xmm1
mulsd   512(%rdi), %xmm1
mulsd   512(%rdi), %xmm1
mulsd   512(%rdi), %xmm1
mulsd   512(%rdi), %xmm1
mulsd   512(%rdi), %xmm1
mulsd   512(%rdi), %xmm1
mulsd   512(%rdi), %xmm1
mulsd   512(%rdi), %xmm1
mulsd   512(%rdi), %xmm1
mulsd   512(%rdi), %xmm1
mulsd   512(%rdi), %xmm1
mulsd   512(%rdi), %xmm1
mulsd   512(%rdi), %xmm1
mulsd   512(%rdi), %xmm1
mulsd   512(%rdi), %xmm1
mulsd   512(%rdi), %xmm1
mulsd   512(%rdi), %xmm1
addsd   %xmm1, %xmm0
ucomisd	%xmm1, %xmm0
setbe	%al
L_end:
rep ret


intersect:
movsd   (%rdi), %xmm0
movsd   %xmm0, (%rdx)
movsd   (%rsi), %xmm0
movsd   %xmm0, (%rcx)

intersect2:

# assume xmm0 has t00, xmm1 has t11
movsd  %xmm0,(%rdx)
movsd  %xmm1,(%rcx)
movabs $0xff00ff00ff00ff00,%rax
vmovq  %rax,%xmm0
vmovq  %rax,%xmm1
vmovq  %rax,%xmm0
vmovq  %rax,%xmm1
vmovq  %rax,%xmm2
vmovq  %rax,%xmm3
vmovq  %rax,%xmm4
vmovq  %rax,%xmm5
vxorpd %xmm2,%xmm2,%xmm2
vxorpd %xmm3,%xmm3,%xmm3
vxorpd %xmm6,%xmm6,%xmm6

vfmadd231sd     8(%rdi), %xmm4, %xmm2
vfmadd231sd     16(%rdi), %xmm4, %xmm2
vfmadd231sd     0x100(%rdi), %xmm4, %xmm2
vfmadd231sd     0x8(%rsi), %xmm4, %xmm3
vfmadd231sd     0x100(%rsi), %xmm4, %xmm3


movaps %xmm2,%xmm0
movaps %xmm3,%xmm1

vsubsd %xmm2, %xmm4, %xmm2 
vdivsd %xmm3, %xmm2, %xmm2 

vcmppd $17, %ymm6, %ymm3, %ymm4 # 17: OP := _CMP_LT_OQ
vcmppd $30, %ymm6, %ymm3, %ymm5 # 30: OP := _CMP_GT_OQ

vcmpsd $17, %xmm6, %xmm3, %xmm4 # 17: OP := _CMP_LT_OQ
vcmpsd $30, %xmm6, %xmm3, %xmm5 # 30: OP := _CMP_GT_OQ

vblendvpd %xmm4, %xmm0, %xmm2, %xmm4
vblendvpd %xmm5, %xmm1, %xmm2, %xmm5

vmaxpd %xmm0, %xmm4, %xmm0
vminpd %xmm1, %xmm5, %xmm1


movaps 0x100(%rsi),%xmm0
movaps (%rsi),%xmm0
movaps %xmm4,%xmm1


movaps %xmm4,%xmm0
movaps %xmm5,%xmm1
movaps %xmm3,%xmm1
movaps %xmm2,%xmm0


vblendvpd %xmm4, %xmm2, %xmm0, %xmm4
vblendvpd %xmm5, %xmm2, %xmm1, %xmm5


vmovq %xmm0,0x100(%rsi)

vmovq  %rax,%xmm0
vmovq  %rax,%xmm1
vmovq  %rax,%xmm2
vmovq  %rax,%xmm3
vmovq  %rax,%xmm4


movsd  %xmm0,(%rsi)
movsd  %xmm1,(%rdx)

nop

lea 0x10(%rip),%rax
mov %edi,%edi
movslq (%rax,%rdi,4),%rcx
add %rax,%rcx
jmpq *%rcx


vxorpd %xmm0,%xmm0,%xmm0
vxorpd %xmm1,%xmm1,%xmm1

vsubsd 0x100(%rcx), %xmm3, %xmm2
vmulpd %xmm2, %xmm4, %xmm2

vmaxpd %xmm0,%xmm2,%xmm0
vminpd %xmm1,%xmm2,%xmm1


movsd  0x100(%rcx), %xmm4

movsd (%rcx),%xmm0
movsd %xmm0,(%rdx)
movsd 8(%rcx),%xmm0
movsd %xmm0,(%rsi)

movsd  %xmm4,(%rsi)
movsd  %xmm3,(%rdx)


movslq (%rax,%rdi,4),%r11
add    %rax,%r11
jmpq   *%r11

movsd %xmm0, (%rsi)

# ------------- three    two    one
# one += three*two
vfmadd231sd     8(%rdi), %xmm4, %xmm2
# one = one*two + three
vfmadd213sd     0x100(%rsi), %xmm0, %xmm4
movsd  %xmm4,0x100(%rsi)


vmovq  %rax,%xmm0
subsd  %xmm1,%xmm0
vmulsd 0x100(%rcx),%xmm4,%xmm2
vfmsub213sd 0x100(%rsi),%xmm0,%xmm4

vperm2f128 $8,%ymm0,%ymm1,%ymm2


vmovq %rax,%xmm0
vmovq %rax,%xmm1
vmovq %rax,%xmm2
vmovq %rax,%xmm3
vmovq %rax,%xmm4
vmovq %rax,%xmm5
vmovq %rax,%xmm6
vmovq %rax,%xmm7
vmovq %rax,%xmm8
vmovq %rax,%xmm9
vmovq %rax,%xmm10
vmovq %rax,%xmm11
vmovq %rax,%xmm12
vmovq %rax,%xmm13
vmovq %rax,%xmm14
vmovq %rax,%xmm15

vmovsd 0x100(%rip), %xmm0
vmovsd 0x100(%rip), %xmm1
vmovsd 0x100(%rip), %xmm2
vmovsd 0x100(%rip), %xmm3
vmovsd 0x100(%rip), %xmm4
vmovsd 0x100(%rip), %xmm5
vmovsd 0x100(%rip), %xmm6
vmovsd 0x100(%rip), %xmm7
vmovsd 0x100(%rip), %xmm8
vmovsd 0x100(%rip), %xmm9
vmovsd 0x100(%rip), %xmm10
vmovsd 0x100(%rip), %xmm11
vmovsd 0x100(%rip), %xmm12
vmovsd 0x100(%rip), %xmm13
vmovsd 0x100(%rip), %xmm14
vmovsd 0x100(%rip), %xmm15

vxorpd %ymm0 ,%ymm0 ,%ymm0
vxorpd %ymm1 ,%ymm1 ,%ymm1
vxorpd %ymm2 ,%ymm2 ,%ymm2
vxorpd %ymm3 ,%ymm3 ,%ymm3
vxorpd %ymm4 ,%ymm4 ,%ymm4
vxorpd %ymm5 ,%ymm5 ,%ymm5
vxorpd %ymm6 ,%ymm6 ,%ymm6
vxorpd %ymm7 ,%ymm7 ,%ymm7
vxorpd %ymm8 ,%ymm8 ,%ymm8
vxorpd %ymm9 ,%ymm9 ,%ymm9
vxorpd %ymm10,%ymm10,%ymm10
vxorpd %ymm11,%ymm11,%ymm11
vxorpd %ymm12,%ymm12,%ymm12
vxorpd %ymm13,%ymm13,%ymm13
vxorpd %ymm14,%ymm14,%ymm14
vxorpd %ymm15,%ymm15,%ymm15


# 1 lat, 1 ipc
# imm, b,a, dst -- 64bit within 128 boundaries
vshufpd $0b10101010,%ymm3,%ymm2,%ymm1

# 1,1
# imm src dst
vpermilpd $0b10101010,%ymm0,%ymm0





