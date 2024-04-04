	.text
	.file	"_vec_not_any_lt"
	.section	.rodata.cst16,"aM",@progbits,16
	.p2align	4                               # -- Begin function julia__vec_not_any_lt_5899
.LCPI0_0:
	.byte	0                               # 0x0
	.byte	1                               # 0x1
	.byte	2                               # 0x2
	.byte	3                               # 0x3
	.byte	4                               # 0x4
	.byte	5                               # 0x5
	.byte	6                               # 0x6
	.byte	7                               # 0x7
	.byte	8                               # 0x8
	.byte	9                               # 0x9
	.byte	10                              # 0xa
	.byte	11                              # 0xb
	.byte	12                              # 0xc
	.byte	13                              # 0xd
	.byte	14                              # 0xe
	.byte	15                              # 0xf
	.text
	.globl	julia__vec_not_any_lt_5899
	.p2align	4, 0x90
	.type	julia__vec_not_any_lt_5899,@function
julia__vec_not_any_lt_5899:             # @julia__vec_not_any_lt_5899
; ┌ @ none within `_vec_not_any_lt`
	.cfi_startproc
# %bb.0:                                # %top
	push	rbp
	.cfi_def_cfa_offset 16
	.cfi_offset rbp, -16
	mov	rbp, rsp
	.cfi_def_cfa_register rbp
; │ @ none within `_vec_not_any_lt` @ none:0
; │┌ @ none within `macro expansion` @ c:\data\projects\gbgb\Groebner.jl\benchmark\scripts\llvm.jl:120
; ││┌ @ abstractarray.jl:1237 within `pointer`
; │││┌ @ pointer.jl:65 within `unsafe_convert`
	mov	r9, qword ptr [rcx]
; ││└└
; ││┌ @ essentials.jl:10 within `length`
	mov	rax, qword ptr [rcx + 8]
; ││└
; ││┌ @ pointer.jl:282 within `+`
	add	r9, 2
; ││└
; ││┌ @ abstractarray.jl:1237 within `pointer`
; │││┌ @ pointer.jl:65 within `unsafe_convert`
	mov	r10, qword ptr [rdx]
; ││└└
; ││┌ @ pointer.jl:282 within `+`
	add	r10, 2
; ││└
; ││┌ @ int.jl:86 within `-`
	lea	r11, [rax - 1]
; ││└
	cmp	r11, 16
	jb	.LBB0_1
# %bb.2:                                # %L9.lr.ph.i
	add	rax, -16
	movabs	r8, 9223372036854775792
	and	r8, r11
	xor	ecx, ecx
	.p2align	4, 0x90
.LBB0_3:                                # %L9.i
                                        # =>This Inner Loop Header: Depth=1
	vmovdqu	ymm0, ymmword ptr [r9 + 2*rcx]
	vpmaxuw	ymm1, ymm0, ymmword ptr [r10 + 2*rcx]
	vpcmpeqw	ymm0, ymm0, ymm1
	vpmovmskb	edx, ymm0
	not	edx
	test	edx, edx
	jne	.LBB0_4
# %bb.5:                                # %L30.i
                                        #   in Loop: Header=BB0_3 Depth=1
	add	rcx, 16
	cmp	rcx, rax
	jl	.LBB0_3
	jmp	.LBB0_6
.LBB0_1:
	xor	r8d, r8d
.LBB0_6:                                # %L32.i
	mov	al, 1
	cmp	r8, r11
	je	.LBB0_72
# %bb.7:                                # %L51.i
	and	r11b, 15
	vmovd	xmm0, r11d
	vpbroadcastb	xmm0, xmm0
	lea	rdx, [r9 + 2*r8]
	movabs	rax, offset .LCPI0_0
	vpcmpgtb	xmm0, xmm0, xmmword ptr [rax]
	vpmovmskb	eax, xmm0
	vpxor	xmm0, xmm0, xmm0
	test	al, 1
	jne	.LBB0_8
# %bb.9:                                # %else
	test	al, 2
	jne	.LBB0_10
.LBB0_11:                               # %else12
	test	al, 4
	jne	.LBB0_12
.LBB0_13:                               # %else15
	test	al, 8
	jne	.LBB0_14
.LBB0_15:                               # %else18
	test	al, 16
	jne	.LBB0_16
.LBB0_17:                               # %else21
	test	al, 32
	jne	.LBB0_18
.LBB0_19:                               # %else24
	test	al, 64
	jne	.LBB0_20
.LBB0_21:                               # %else27
	test	al, -128
	jne	.LBB0_22
.LBB0_23:                               # %else30
	test	eax, 256
	jne	.LBB0_24
.LBB0_25:                               # %else33
	test	eax, 512
	jne	.LBB0_26
.LBB0_27:                               # %else36
	test	eax, 1024
	jne	.LBB0_28
.LBB0_29:                               # %else39
	test	eax, 2048
	jne	.LBB0_30
.LBB0_31:                               # %else42
	test	eax, 4096
	jne	.LBB0_32
.LBB0_33:                               # %else45
	test	eax, 8192
	jne	.LBB0_34
.LBB0_35:                               # %else48
	test	eax, 16384
	jne	.LBB0_36
.LBB0_37:                               # %else51
	test	eax, 32768
	je	.LBB0_39
.LBB0_38:                               # %cond.load53
	vpbroadcastw	ymm1, word ptr [rdx + 30]
	vpblendw	ymm1, ymm0, ymm1, 128           # ymm1 = ymm0[0,1,2,3,4,5,6],ymm1[7],ymm0[8,9,10,11,12,13,14],ymm1[15]
	vpblendd	ymm0, ymm0, ymm1, 240           # ymm0 = ymm0[0,1,2,3],ymm1[4,5,6,7]
.LBB0_39:                               # %else54
	lea	rcx, [r10 + 2*r8]
	vpxor	xmm1, xmm1, xmm1
	test	al, 1
	jne	.LBB0_40
# %bb.41:                               # %else58
	test	al, 2
	jne	.LBB0_42
.LBB0_43:                               # %else61
	test	al, 4
	jne	.LBB0_44
.LBB0_45:                               # %else64
	test	al, 8
	jne	.LBB0_46
.LBB0_47:                               # %else67
	test	al, 16
	jne	.LBB0_48
.LBB0_49:                               # %else70
	test	al, 32
	jne	.LBB0_50
.LBB0_51:                               # %else73
	test	al, 64
	jne	.LBB0_52
.LBB0_53:                               # %else76
	test	al, -128
	jne	.LBB0_54
.LBB0_55:                               # %else79
	test	eax, 256
	jne	.LBB0_56
.LBB0_57:                               # %else82
	test	eax, 512
	jne	.LBB0_58
.LBB0_59:                               # %else85
	test	eax, 1024
	jne	.LBB0_60
.LBB0_61:                               # %else88
	test	eax, 2048
	jne	.LBB0_62
.LBB0_63:                               # %else91
	test	eax, 4096
	jne	.LBB0_64
.LBB0_65:                               # %else94
	test	eax, 8192
	jne	.LBB0_66
.LBB0_67:                               # %else97
	test	eax, 16384
	jne	.LBB0_68
.LBB0_69:                               # %else100
	test	eax, 32768
	je	.LBB0_71
.LBB0_70:                               # %cond.load102
	vpbroadcastw	ymm2, word ptr [rcx + 30]
	vpblendw	ymm2, ymm1, ymm2, 128           # ymm2 = ymm1[0,1,2,3,4,5,6],ymm2[7],ymm1[8,9,10,11,12,13,14],ymm2[15]
	vpblendd	ymm1, ymm1, ymm2, 240           # ymm1 = ymm1[0,1,2,3],ymm2[4,5,6,7]
.LBB0_71:                               # %else103
	vpmaxuw	ymm1, ymm0, ymm1
	vpcmpeqw	ymm0, ymm0, ymm1
	vpmovmskb	eax, ymm0
	not	eax
	test	eax, eax
	sete	al
.LBB0_72:                               # %julia__vec_not_any_lt_5899u5901.exit
                                        # kill: def $al killed $al killed $eax
	pop	rbp
	vzeroupper
	ret
.LBB0_4:
	xor	eax, eax
                                        # kill: def $al killed $al killed $eax
	pop	rbp
	vzeroupper
	ret
.LBB0_8:                                # %cond.load
	movzx	ecx, word ptr [rdx]
	vmovd	xmm0, ecx
	test	al, 2
	je	.LBB0_11
.LBB0_10:                               # %cond.load11
	vpinsrw	xmm1, xmm0, word ptr [rdx + 2], 1
	vpblendd	ymm0, ymm0, ymm1, 15            # ymm0 = ymm1[0,1,2,3],ymm0[4,5,6,7]
	test	al, 4
	je	.LBB0_13
.LBB0_12:                               # %cond.load14
	vpinsrw	xmm1, xmm0, word ptr [rdx + 4], 2
	vpblendd	ymm0, ymm0, ymm1, 15            # ymm0 = ymm1[0,1,2,3],ymm0[4,5,6,7]
	test	al, 8
	je	.LBB0_15
.LBB0_14:                               # %cond.load17
	vpinsrw	xmm1, xmm0, word ptr [rdx + 6], 3
	vpblendd	ymm0, ymm0, ymm1, 15            # ymm0 = ymm1[0,1,2,3],ymm0[4,5,6,7]
	test	al, 16
	je	.LBB0_17
.LBB0_16:                               # %cond.load20
	vpinsrw	xmm1, xmm0, word ptr [rdx + 8], 4
	vpblendd	ymm0, ymm0, ymm1, 15            # ymm0 = ymm1[0,1,2,3],ymm0[4,5,6,7]
	test	al, 32
	je	.LBB0_19
.LBB0_18:                               # %cond.load23
	vpinsrw	xmm1, xmm0, word ptr [rdx + 10], 5
	vpblendd	ymm0, ymm0, ymm1, 15            # ymm0 = ymm1[0,1,2,3],ymm0[4,5,6,7]
	test	al, 64
	je	.LBB0_21
.LBB0_20:                               # %cond.load26
	vpinsrw	xmm1, xmm0, word ptr [rdx + 12], 6
	vpblendd	ymm0, ymm0, ymm1, 15            # ymm0 = ymm1[0,1,2,3],ymm0[4,5,6,7]
	test	al, -128
	je	.LBB0_23
.LBB0_22:                               # %cond.load29
	vpinsrw	xmm1, xmm0, word ptr [rdx + 14], 7
	vpblendd	ymm0, ymm0, ymm1, 15            # ymm0 = ymm1[0,1,2,3],ymm0[4,5,6,7]
	test	eax, 256
	je	.LBB0_25
.LBB0_24:                               # %cond.load32
	vpbroadcastw	ymm1, word ptr [rdx + 16]
	vpblendw	ymm1, ymm0, ymm1, 1             # ymm1 = ymm1[0],ymm0[1,2,3,4,5,6,7],ymm1[8],ymm0[9,10,11,12,13,14,15]
	vpblendd	ymm0, ymm0, ymm1, 240           # ymm0 = ymm0[0,1,2,3],ymm1[4,5,6,7]
	test	eax, 512
	je	.LBB0_27
.LBB0_26:                               # %cond.load35
	vpbroadcastw	ymm1, word ptr [rdx + 18]
	vpblendw	ymm1, ymm0, ymm1, 2             # ymm1 = ymm0[0],ymm1[1],ymm0[2,3,4,5,6,7,8],ymm1[9],ymm0[10,11,12,13,14,15]
	vpblendd	ymm0, ymm0, ymm1, 240           # ymm0 = ymm0[0,1,2,3],ymm1[4,5,6,7]
	test	eax, 1024
	je	.LBB0_29
.LBB0_28:                               # %cond.load38
	vpbroadcastw	ymm1, word ptr [rdx + 20]
	vpblendw	ymm1, ymm0, ymm1, 4             # ymm1 = ymm0[0,1],ymm1[2],ymm0[3,4,5,6,7,8,9],ymm1[10],ymm0[11,12,13,14,15]
	vpblendd	ymm0, ymm0, ymm1, 240           # ymm0 = ymm0[0,1,2,3],ymm1[4,5,6,7]
	test	eax, 2048
	je	.LBB0_31
.LBB0_30:                               # %cond.load41
	vpbroadcastw	ymm1, word ptr [rdx + 22]
	vpblendw	ymm1, ymm0, ymm1, 8             # ymm1 = ymm0[0,1,2],ymm1[3],ymm0[4,5,6,7,8,9,10],ymm1[11],ymm0[12,13,14,15]
	vpblendd	ymm0, ymm0, ymm1, 240           # ymm0 = ymm0[0,1,2,3],ymm1[4,5,6,7]
	test	eax, 4096
	je	.LBB0_33
.LBB0_32:                               # %cond.load44
	vpbroadcastw	ymm1, word ptr [rdx + 24]
	vpblendw	ymm1, ymm0, ymm1, 16            # ymm1 = ymm0[0,1,2,3],ymm1[4],ymm0[5,6,7,8,9,10,11],ymm1[12],ymm0[13,14,15]
	vpblendd	ymm0, ymm0, ymm1, 240           # ymm0 = ymm0[0,1,2,3],ymm1[4,5,6,7]
	test	eax, 8192
	je	.LBB0_35
.LBB0_34:                               # %cond.load47
	vpbroadcastw	ymm1, word ptr [rdx + 26]
	vpblendw	ymm1, ymm0, ymm1, 32            # ymm1 = ymm0[0,1,2,3,4],ymm1[5],ymm0[6,7,8,9,10,11,12],ymm1[13],ymm0[14,15]
	vpblendd	ymm0, ymm0, ymm1, 240           # ymm0 = ymm0[0,1,2,3],ymm1[4,5,6,7]
	test	eax, 16384
	je	.LBB0_37
.LBB0_36:                               # %cond.load50
	vpbroadcastw	ymm1, word ptr [rdx + 28]
	vpblendw	ymm1, ymm0, ymm1, 64            # ymm1 = ymm0[0,1,2,3,4,5],ymm1[6],ymm0[7,8,9,10,11,12,13],ymm1[14],ymm0[15]
	vpblendd	ymm0, ymm0, ymm1, 240           # ymm0 = ymm0[0,1,2,3],ymm1[4,5,6,7]
	test	eax, 32768
	je	.LBB0_39
	jmp	.LBB0_38
.LBB0_40:                               # %cond.load57
	movzx	edx, word ptr [rcx]
	vmovd	xmm1, edx
	test	al, 2
	je	.LBB0_43
.LBB0_42:                               # %cond.load60
	vpinsrw	xmm2, xmm1, word ptr [rcx + 2], 1
	vpblendd	ymm1, ymm1, ymm2, 15            # ymm1 = ymm2[0,1,2,3],ymm1[4,5,6,7]
	test	al, 4
	je	.LBB0_45
.LBB0_44:                               # %cond.load63
	vpinsrw	xmm2, xmm1, word ptr [rcx + 4], 2
	vpblendd	ymm1, ymm1, ymm2, 15            # ymm1 = ymm2[0,1,2,3],ymm1[4,5,6,7]
	test	al, 8
	je	.LBB0_47
.LBB0_46:                               # %cond.load66
	vpinsrw	xmm2, xmm1, word ptr [rcx + 6], 3
	vpblendd	ymm1, ymm1, ymm2, 15            # ymm1 = ymm2[0,1,2,3],ymm1[4,5,6,7]
	test	al, 16
	je	.LBB0_49
.LBB0_48:                               # %cond.load69
	vpinsrw	xmm2, xmm1, word ptr [rcx + 8], 4
	vpblendd	ymm1, ymm1, ymm2, 15            # ymm1 = ymm2[0,1,2,3],ymm1[4,5,6,7]
	test	al, 32
	je	.LBB0_51
.LBB0_50:                               # %cond.load72
	vpinsrw	xmm2, xmm1, word ptr [rcx + 10], 5
	vpblendd	ymm1, ymm1, ymm2, 15            # ymm1 = ymm2[0,1,2,3],ymm1[4,5,6,7]
	test	al, 64
	je	.LBB0_53
.LBB0_52:                               # %cond.load75
	vpinsrw	xmm2, xmm1, word ptr [rcx + 12], 6
	vpblendd	ymm1, ymm1, ymm2, 15            # ymm1 = ymm2[0,1,2,3],ymm1[4,5,6,7]
	test	al, -128
	je	.LBB0_55
.LBB0_54:                               # %cond.load78
	vpinsrw	xmm2, xmm1, word ptr [rcx + 14], 7
	vpblendd	ymm1, ymm1, ymm2, 15            # ymm1 = ymm2[0,1,2,3],ymm1[4,5,6,7]
	test	eax, 256
	je	.LBB0_57
.LBB0_56:                               # %cond.load81
	vpbroadcastw	ymm2, word ptr [rcx + 16]
	vpblendw	ymm2, ymm1, ymm2, 1             # ymm2 = ymm2[0],ymm1[1,2,3,4,5,6,7],ymm2[8],ymm1[9,10,11,12,13,14,15]
	vpblendd	ymm1, ymm1, ymm2, 240           # ymm1 = ymm1[0,1,2,3],ymm2[4,5,6,7]
	test	eax, 512
	je	.LBB0_59
.LBB0_58:                               # %cond.load84
	vpbroadcastw	ymm2, word ptr [rcx + 18]
	vpblendw	ymm2, ymm1, ymm2, 2             # ymm2 = ymm1[0],ymm2[1],ymm1[2,3,4,5,6,7,8],ymm2[9],ymm1[10,11,12,13,14,15]
	vpblendd	ymm1, ymm1, ymm2, 240           # ymm1 = ymm1[0,1,2,3],ymm2[4,5,6,7]
	test	eax, 1024
	je	.LBB0_61
.LBB0_60:                               # %cond.load87
	vpbroadcastw	ymm2, word ptr [rcx + 20]
	vpblendw	ymm2, ymm1, ymm2, 4             # ymm2 = ymm1[0,1],ymm2[2],ymm1[3,4,5,6,7,8,9],ymm2[10],ymm1[11,12,13,14,15]
	vpblendd	ymm1, ymm1, ymm2, 240           # ymm1 = ymm1[0,1,2,3],ymm2[4,5,6,7]
	test	eax, 2048
	je	.LBB0_63
.LBB0_62:                               # %cond.load90
	vpbroadcastw	ymm2, word ptr [rcx + 22]
	vpblendw	ymm2, ymm1, ymm2, 8             # ymm2 = ymm1[0,1,2],ymm2[3],ymm1[4,5,6,7,8,9,10],ymm2[11],ymm1[12,13,14,15]
	vpblendd	ymm1, ymm1, ymm2, 240           # ymm1 = ymm1[0,1,2,3],ymm2[4,5,6,7]
	test	eax, 4096
	je	.LBB0_65
.LBB0_64:                               # %cond.load93
	vpbroadcastw	ymm2, word ptr [rcx + 24]
	vpblendw	ymm2, ymm1, ymm2, 16            # ymm2 = ymm1[0,1,2,3],ymm2[4],ymm1[5,6,7,8,9,10,11],ymm2[12],ymm1[13,14,15]
	vpblendd	ymm1, ymm1, ymm2, 240           # ymm1 = ymm1[0,1,2,3],ymm2[4,5,6,7]
	test	eax, 8192
	je	.LBB0_67
.LBB0_66:                               # %cond.load96
	vpbroadcastw	ymm2, word ptr [rcx + 26]
	vpblendw	ymm2, ymm1, ymm2, 32            # ymm2 = ymm1[0,1,2,3,4],ymm2[5],ymm1[6,7,8,9,10,11,12],ymm2[13],ymm1[14,15]
	vpblendd	ymm1, ymm1, ymm2, 240           # ymm1 = ymm1[0,1,2,3],ymm2[4,5,6,7]
	test	eax, 16384
	je	.LBB0_69
.LBB0_68:                               # %cond.load99
	vpbroadcastw	ymm2, word ptr [rcx + 28]
	vpblendw	ymm2, ymm1, ymm2, 64            # ymm2 = ymm1[0,1,2,3,4,5],ymm2[6],ymm1[7,8,9,10,11,12,13],ymm2[14],ymm1[15]
	vpblendd	ymm1, ymm1, ymm2, 240           # ymm1 = ymm1[0,1,2,3],ymm2[4,5,6,7]
	test	eax, 32768
	je	.LBB0_71
	jmp	.LBB0_70
.Lfunc_end0:
	.size	julia__vec_not_any_lt_5899, .Lfunc_end0-julia__vec_not_any_lt_5899
	.cfi_endproc
; └└
                                        # -- End function
	.section	".note.GNU-stack","",@progbits
