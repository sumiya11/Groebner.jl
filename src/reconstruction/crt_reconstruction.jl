# CRT

"""
    CRT!

Implements the Chinese Remainder lifting algorithm. Computes the unique `x` such
that
        
    x ≡ a1 mod m1
    x ≡ a2 mod m2

Writes the answer to `buf` inplace.

## Additional parameters:

- `M`: must be equal to `m1 * m2`
- `buf`, `n1`, `n2`: additional buffers
- `c1`, `c2`: numbers, such that `c1 = m2 * invmod(m2, m1)`, and `c2 = m1 *
  invmod(m1, m2)`

Then, `x` is obtained as `x = c1 a1 + c2 a2 mod M`.
"""
function CRT!(
    M::BigInt,
    buf::BigInt,
    n1::BigInt,
    n2::BigInt,
    a1::BigInt,
    c1::BigInt,
    a2::UInt,
    c2::BigInt
)
    Base.GMP.MPZ.mul_ui!(n1, c2, a2)
    Base.GMP.MPZ.mul!(n2, c1, a1)

    Base.GMP.MPZ.add!(buf, n1, n2)
    Base.GMP.MPZ.fdiv_r!(buf, M)

    nothing
end

"""
    CRT_precompute!

Given the two moduli `m1` and `m2`, precomputes the multipliers `c1`, `c2` and
the modulo `M` for CRT.

`n1` and `n2` are additional buffers.
"""
function CRT_precompute!(
    M::BigInt,
    n1::BigInt,
    n2::BigInt,
    c1::BigInt,
    m1::BigInt,
    c2::BigInt,
    m2::UInt
)
    Base.GMP.MPZ.mul_ui!(M, m1, m2)
    Base.GMP.MPZ.set_ui!(n2, m2)

    Base.GMP.MPZ.gcdext!(n1, c1, c2, n2, m1)
    Base.GMP.MPZ.mul_ui!(c1, m2)
    Base.GMP.MPZ.mul!(c2, m1)

    nothing
end

function CRT!(
    M::BigInt,
    buf::BigInt,
    n1::BigInt,
    n2::BigInt,
    a1::BigInt,
    c1::BigInt,
    a2::BigInt,
    c2::BigInt
)
    Base.GMP.MPZ.mul!(n1, c2, a2)
    Base.GMP.MPZ.mul!(n2, c1, a1)

    Base.GMP.MPZ.add!(buf, n1, n2)
    Base.GMP.MPZ.fdiv_r!(buf, M)

    nothing
end

function CRT_precompute!(
    M::BigInt,
    n1::BigInt,
    n2::BigInt,
    c1::BigInt,
    m1::BigInt,
    c2::BigInt,
    m2::BigInt
)
    Base.GMP.MPZ.mul!(M, m1, m2)

    Base.GMP.MPZ.gcdext!(n1, c1, c2, m2, m1)
    Base.GMP.MPZ.mul!(c1, m2)
    Base.GMP.MPZ.mul!(c2, m1)

    nothing
end

# CRT for multiple remainders
function CRT!(
    M::BigInt,
    buf::BigInt,
    n1::BigInt,
    n2::BigInt,
    ai::Vector{UInt},
    ci::Vector{BigInt}
)
    @invariant length(ai) == length(ci)

    Base.GMP.MPZ.set_ui!(n1, UInt(0))
    for i in 1:length(ai)
        Base.GMP.MPZ.mul_ui!(n2, ci[i], ai[i])
        Base.GMP.MPZ.add!(n1, n2)
    end

    Base.GMP.MPZ.set!(buf, n1)
    Base.GMP.MPZ.fdiv_r!(buf, M)

    nothing
end

function CRT_precompute!(
    M::BigInt,
    n1::BigInt,
    n2::BigInt,
    ci::Vector{BigInt},
    moduli::Vector{UInt}
)
    @invariant length(ci) == length(moduli)

    n3, n4 = BigInt(), BigInt()
    Base.GMP.MPZ.set_ui!(M, moduli[1])
    @inbounds for i in 2:length(moduli)
        Base.GMP.MPZ.mul_ui!(M, moduli[i])
    end

    @inbounds for i in 1:length(moduli)
        Base.GMP.MPZ.set_ui!(n2, moduli[i])
        Base.GMP.MPZ.tdiv_q!(ci[i], M, n2)
        Base.GMP.MPZ.gcdext!(n2, n3, n4, ci[i], n2)
        Base.GMP.MPZ.mul!(ci[i], n3)
    end

    nothing
end
