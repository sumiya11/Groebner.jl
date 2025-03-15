# This file is a part of Groebner.jl. License is GNU GPL v2.

###
# CRT

# Auxiliary functions allowing to avoid overflow
# on Windows for numbers between 32 and 64 bits
function my_set_ui!(a::BigInt, b::UInt64)
    @static if Culong == UInt64
        Base.GMP.MPZ.set_ui!(a, b)
    else
        if b < 2^32
            Base.GMP.MPZ.set_ui!(a, b)
        else
            Base.GMP.MPZ.set!(a, BigInt(b))
        end
    end
end

function my_mul_ui!(a::BigInt, b::BigInt, c::UInt64)
    @static if Culong == UInt64
        Base.GMP.MPZ.mul_ui!(a, b, c)
    else
        if c < 2^32
            Base.GMP.MPZ.mul_ui!(a, b, c)
        else
            Base.GMP.MPZ.mul!(a, b, BigInt(c))
        end
    end
end

function my_mul_ui!(a::BigInt, b::UInt64)
    @static if Culong == UInt64
        Base.GMP.MPZ.mul_ui!(a, b)
    else
        if b < 2^32
            Base.GMP.MPZ.mul_ui!(a, b)
        else
            Base.GMP.MPZ.mul!(a, BigInt(b))
        end
    end
end

"""
    crt!

Chinese Remainder lifting. Computes the unique `x` such that for all `i`
        
    x ≡ ai mod mi.

Writes the answer to `buf` inplace.

## Additional parameters:

- `M`: must be equal to `∏ mi`
- `buf`, `n1`, `n2`: additional buffers
- `ci`: an array of numbers, with `ci[i] = πi * invmod(πi, mi)`, `πi = M / mi`

Then, `x` is obtained as `x = ∑ ci[i] ai[i] mod M`.
"""
function crt!(M::BigInt, buf::BigInt, n1::BigInt, n2::BigInt, ai::Vector{UInt}, ci::Vector{BigInt})
    @invariant length(ai) == length(ci)

    my_set_ui!(n1, UInt(0))
    for i in 1:length(ai)
        my_mul_ui!(n2, ci[i], ai[i])
        Base.GMP.MPZ.add!(n1, n2)
    end

    Base.GMP.MPZ.set!(buf, n1)
    Base.GMP.MPZ.fdiv_r!(buf, M)

    nothing
end

"""
    crt_precompute!

Given the moduli, precomputes the multipliers `ci`, and the modulo `M` for CRT.

`n1` and `n2` are additional buffers.
"""
function crt_precompute!(
    M::BigInt,
    n1::BigInt,
    n2::BigInt,
    ci::Vector{BigInt},
    moduli::Vector{UInt}
)
    @invariant length(ci) == length(moduli)

    n3, n4 = BigInt(), BigInt()
    @inbounds my_set_ui!(M, moduli[1])
    @inbounds for i in 2:length(moduli)
        my_mul_ui!(M, moduli[i])
    end

    @inbounds for i in 1:length(moduli)
        my_set_ui!(n2, moduli[i])
        Base.GMP.MPZ.tdiv_q!(ci[i], M, n2)
        Base.GMP.MPZ.gcdext!(n2, n3, n4, ci[i], n2)
        Base.GMP.MPZ.mul!(ci[i], n3)
    end

    nothing
end

###
# Element-wise CRT

# Reconstructs only the witness set.
# Table of big integers must be initialized.
function crt_vec_partial!(
    table_zz::Vector{Vector{BigInt}},
    modulo::BigInt,
    tables_ff::Vector{Vector{Vector{T}}},
    moduli::Vector{U},
    witness_set::Vector{Tuple{Int, Int}},
    mask::Vector{BitVector}
) where {T <: Integer, U <: Integer}
    @invariant length(tables_ff) == length(moduli)
    @invariant length(table_zz) == length(mask)

    # Precompute CRT multipliers
    buf, n1, n2 = BigInt(), BigInt(), BigInt()
    mults = Vector{BigInt}(undef, length(moduli))
    for i in 1:length(mults)
        mults[i] = BigInt(0)
    end

    crt_precompute!(modulo, n1, n2, mults, map(UInt64, moduli))

    rems = Vector{UInt64}(undef, length(moduli))
    @inbounds for k in 1:length(witness_set)
        i, j = witness_set[k]

        for t in 1:length(moduli)
            rems[t] = UInt64(tables_ff[t][i][j])
        end

        crt!(modulo, buf, n1, n2, rems, mults)

        Base.GMP.MPZ.set!(table_zz[i][j], buf)

        mask[i][j] = true
    end

    nothing
end

# Table of big integers must be initialized.
function crt_vec_full!(
    table_zz::Vector{Vector{BigInt}},
    modulo::BigInt,
    tables_ff::Vector{Vector{Vector{T}}},
    moduli::Vector{U},
    mask::Vector{BitVector}
) where {T <: Integer, U <: Integer}
    @invariant length(tables_ff) == length(moduli)
    @invariant length(table_zz) == length(mask)

    buf, n1, n2 = BigInt(), BigInt(), BigInt()
    mults = Vector{BigInt}(undef, length(moduli))
    for i in 1:length(mults)
        mults[i] = BigInt(0)
    end

    crt_precompute!(modulo, n1, n2, mults, map(UInt64, moduli))

    rems = Vector{UInt64}(undef, length(moduli))
    @inbounds for i in 1:length(table_zz)
        for j in 1:length(table_zz[i])
            mask[i][j] && continue

            for k in 1:length(moduli)
                @invariant 0 <= tables_ff[k][i][j] < moduli[k]
                rems[k] = UInt64(tables_ff[k][i][j])
            end

            crt!(modulo, buf, n1, n2, rems, mults)

            Base.GMP.MPZ.set!(table_zz[i][j], buf)
        end
    end
end
