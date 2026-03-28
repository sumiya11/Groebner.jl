## Benchmarking monomial implementations

### Generating native code

```julia
julia native.jl
```

### Benchmarking

```julia
julia benchmark.jl
```

### Findings

i5-8250U CPU @ 1.60GHz (2017):

1. `monom_hash` is not inlined for `PackedTuple2`. Also, `monom_hash` checks inbounds in `view` in `PackedTuple2`.

2. `monom_create_divmask` is inlined for `FixedMonomNoDeg` (and not inlined for `Vector{UInt8}` or `PackedTuple2`). Produces +400 lines of assembly => **mark @noinline?**.

3. As a result of 2., xmm6-xmm15 are spilled per calling convention. If `hashtable_insert` is not inlined itself, then the spills may be expensive (=> **mark `monom_create_divmask` @noinline?**).

4. A notable difference between implementations is whether the total degree contributes to hash.

5. Combining exponents in `monom_hash` of `NibbleMonom` is fine.

6. Inserting in the hashtable (`hashtable_insert!`) many times:

- Original code (with `monom_create_divmask` inlined)
```
=======================================================
 Benchmark Results for N = 11 (Min ns / Insertion)
=======================================================
Monomial Type                  | Empty HT   | Filled HT
-------------------------------+------------+----------
PackedTuple2{UInt64, UInt8}    |   104.29   |    22.00
Vector{UInt8}                  |   118.66   |    26.59
FixedVector{16, UInt8}         |    78.90   |    16.51
FixedMonom{16, UInt8}          |    83.20   |    18.08
FixedMonomNoDeg{16, UInt8}     |    81.03   |    16.11
NibbleMonom{8}                 |    82.31   |    11.26
=======================================================

=======================================================
 Benchmark Results for N = 55 (Min ns / Insertion)
=======================================================
Monomial Type                  | Empty HT   | Filled HT
-------------------------------+------------+----------
Vector{UInt8}                  |   164.00   |    41.89
FixedVector{64, UInt8}         |   157.22   |    25.47
FixedMonom{64, UInt8}          |   177.21   |    34.78
FixedMonomNoDeg{64, UInt8}     |   157.84   |    26.58
NibbleMonom{32}                |   158.18   |    18.68
=======================================================
```

- With `monom_create_divmask` marked @noinline:
```
=======================================================
 Benchmark Results for N = 11 (Min ns / Insertion)
=======================================================
Monomial Type                  | Empty HT   | Filled HT
-------------------------------+------------+----------
PackedTuple2{UInt64, UInt8}    |    90.16   |    21.30
Vector{UInt8}                  |   102.22   |    24.55
FixedVector{16, UInt8}         |    75.57   |    12.55
FixedMonom{16, UInt8}          |    76.25   |    13.09
FixedMonomNoDeg{16, UInt8}     |    75.86   |    12.49
NibbleMonom{8}                 |    69.18   |    11.30
=======================================================

=======================================================
 Benchmark Results for N = 55 (Min ns / Insertion)
=======================================================
Monomial Type                  | Empty HT   | Filled HT
-------------------------------+------------+----------
Vector{UInt8}                  |   148.70   |    38.96
FixedVector{64, UInt8}         |   140.58   |    23.31
FixedMonom{64, UInt8}          |   142.83   |    25.57
FixedMonomNoDeg{64, UInt8}     |   140.79   |    22.75
NibbleMonom{32}                |   129.29   |    19.14
=======================================================
```

### Other

1. `monom_product!` for `FixedMonomNoDeg` for `N=11` (rounded to 16) is nice (checked arithmetic):

```julia
julia> code_native(Groebner.monom_product!, Tuple{Groebner.FixedMonomNoDeg{16,UInt8},Groebner.FixedMonomNoDeg{16,UInt8},Groebner.FixedMonomNoDeg{16,UInt8}}; debuginfo=:none)

        .text
        .file   "monom_product!"
        .section        .ltext,"axl",@progbits
        .globl  "julia_monom_product!_278751"   # -- Begin function julia_monom_product!_278751
        .p2align        4, 0x90
        .type   "julia_monom_product!_278751",@function
"julia_monom_product!_278751":          # @"julia_monom_product!_278751"
; Function Signature: monom_product!(Groebner.FixedMonomNoDeg{16, UInt8}, Groebner.FixedMonomNoDeg{16, UInt8}, Groebner.FixedMonomNoDeg{16, UInt8})
        .cfi_startproc
# %bb.0:                                # %top
        #DEBUG_VALUE: monom_product!:a <- [$r8+0]
        #DEBUG_VALUE: monom_product!:b <- [$r9+0]
        vmovdqu xmm1, xmmword ptr [r8]
        vpaddb  xmm0, xmm1, xmmword ptr [r9]
        vpmaxub xmm1, xmm0, xmm1
        vpxor   xmm1, xmm0, xmm1
        vptest  xmm1, xmm1
        jne     .LBB0_2
# %bb.1:                                # %L193
        vmovdqu xmmword ptr [rcx], xmm0
        mov     rax, rcx
        ret
.LBB0_2:                                # %L191
        push    rbp
        .cfi_def_cfa_offset 16
        .cfi_offset rbp, -16
        mov     rbp, rsp
        .cfi_def_cfa_register rbp
        sub     rsp, 32
        movabs  rcx, offset ".Ljl_global#278758.jit"
        movabs  rax, offset j_overflowerror_278757
        call    rax
.Lfunc_end0:
        .size   "julia_monom_product!_278751", .Lfunc_end0-"julia_monom_product!_278751"
        .cfi_endproc
                                        # -- End function
.set ".L+Groebner.FixedMonomNoDeg#278753.jit", 1806300472528
        .size   ".L+Groebner.FixedMonomNoDeg#278753.jit", 8
.set ".Ljl_global#278758.jit", 140735697973696
        .size   ".Ljl_global#278758.jit", 8
        .section        ".note.GNU-stack","",@progbits
```
