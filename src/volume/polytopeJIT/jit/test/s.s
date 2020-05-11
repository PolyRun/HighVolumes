.global _inside

.text
_inside:
xorl %eax, %eax
vmovq  %rsi,%xmm0
vmovq  %rsi,%xmm1
vmovq  %rsi,%xmm2
ja L_end
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












