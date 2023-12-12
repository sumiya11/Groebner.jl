using InteractiveUtils, BenchmarkTools

function add_mul!(vres::Vector{Int64}, a::Int32, v::Vector{Int32})
    @inbounds for i in eachindex(vres)
        vres[i] += Int64(a) * Int64(v[i])
    end
end

function add_mul_f!(vres, a, v)
    @fastmath @inbounds @simd for i in eachindex(vres)
        vres[i] += a * v[i]
    end
end

n = 2^10
@benchmark add_mul!(v1, c, v2) setup = begin
    v1 = rand(Int, n)
    c = rand(Int32)
    v2 = rand(Int32, n)
end

@benchmark add_mul_f!(v1, c, v2) setup = begin
    v1 = rand(Float64, n)
    c = rand(Float64)
    v2 = rand(Float64, n)
end

@code_native debuginfo = :none add_mul_f!(Float64[], Float32(1), Float32[])

@code_native debuginfo = :none add_mul!(Int64[], Int32(1), Int32[])
#=
        .text
        .file   "add_mul"
        .globl  julia_add_mul_13693             # -- Begin function julia_add_mul_13693
        .p2align        4, 0x90
        .type   julia_add_mul_13693,@function
julia_add_mul_13693:                    # @julia_add_mul_13693
        .cfi_startproc
# %bb.0:                                # %top
        pushq   %rbp
        .cfi_def_cfa_offset 16
        .cfi_offset %rbp, -16
        movq    %rsp, %rbp
        .cfi_def_cfa_register %rbp
        movq    8(%rcx), %r10
        testq   %r10, %r10
        je      .LBB0_9
# %bb.1:                                # %L18.preheader
        movq    (%rcx), %rcx
        movslq  %edx, %r9
        movq    (%r8), %r11
        movl    $1, %eax
        cmpq    $16, %r10
        jb      .LBB0_7
# %bb.2:                                # %vector.memcheck
        leaq    (%r11,%r10,4), %rdx
        cmpq    %rdx, %rcx
        jae     .LBB0_4
# %bb.3:                                # %vector.memcheck
        leaq    (%rcx,%r10,8), %rdx
        cmpq    %rdx, %r11
        jb      .LBB0_7
.LBB0_4:                                # %vector.ph
        movq    %r10, %r8
        andq    $-16, %r8
        leaq    1(%r8), %rax
        vmovq   %r9, %xmm0
        vpbroadcastq    %xmm0, %ymm0
        xorl    %edx, %edx
        .p2align        4, 0x90
.LBB0_5:                                # %vector.body
                                        # =>This Inner Loop Header: Depth=1
        vpmovzxdq       (%r11,%rdx,4), %ymm1    # ymm1 = mem[0],zero,mem[1],zero,mem[2],zero,mem[3],zero
        vpmovzxdq       16(%r11,%rdx,4), %ymm2  # ymm2 = mem[0],zero,mem[1],zero,mem[2],zero,mem[3],zero
        vpmovzxdq       32(%r11,%rdx,4), %ymm3  # ymm3 = mem[0],zero,mem[1],zero,mem[2],zero,mem[3],zero
        vpmovzxdq       48(%r11,%rdx,4), %ymm4  # ymm4 = mem[0],zero,mem[1],zero,mem[2],zero,mem[3],zero
        vpmuldq %ymm1, %ymm0, %ymm1
        vpmuldq %ymm2, %ymm0, %ymm2
        vpmuldq %ymm3, %ymm0, %ymm3
        vpmuldq %ymm4, %ymm0, %ymm4
        vpaddq  (%rcx,%rdx,8), %ymm1, %ymm1
        vpaddq  32(%rcx,%rdx,8), %ymm2, %ymm2
        vpaddq  64(%rcx,%rdx,8), %ymm3, %ymm3
        vpaddq  96(%rcx,%rdx,8), %ymm4, %ymm4
        vmovdqu %ymm1, (%rcx,%rdx,8)
        vmovdqu %ymm2, 32(%rcx,%rdx,8)
        vmovdqu %ymm3, 64(%rcx,%rdx,8)
        vmovdqu %ymm4, 96(%rcx,%rdx,8)
        addq    $16, %rdx
        cmpq    %rdx, %r8
        jne     .LBB0_5
# %bb.6:                                # %middle.block
        cmpq    %r8, %r10
        je      .LBB0_9
.LBB0_7:                                # %scalar.ph
        decq    %rax
        .p2align        4, 0x90
.LBB0_8:                                # %L18
                                        # =>This Inner Loop Header: Depth=1
        movslq  (%r11,%rax,4), %rdx
        imulq   %r9, %rdx
        addq    %rdx, (%rcx,%rax,8)
        incq    %rax
        cmpq    %rax, %r10
        jne     .LBB0_8
.LBB0_9:                                # %L38
        popq    %rbp
        vzeroupper
        retq
.Lfunc_end0:
        .size   julia_add_mul_13693, .Lfunc_end0-julia_add_mul_13693     
        .cfi_endproc
                                        # -- End function
        .section        ".note.GNU-stack","",@progbits
=#
