# Rational reconstruction and Chinese remainder theorem (CRT) reconstruction

"""
Returns the bound for rational reconstruction (see `rational_reconstruction!`)
based on the bitsize of the modulo. As soon as the numerator in rational
reconstruction exceeds this bound, the gcd iteration stops
"""
function rational_reconstruction_bound(modulo::BigInt)
    setprecision(2 * Base.GMP.MPZ.sizeinbase(modulo, 2)) do
        ceil(BigInt, sqrt(BigFloat(modulo) / 2))
    end
end

"""
Computes the rational reconstruction of `a` mod `m`. Namely, a pair of numbers
`num`, `den`, such that 

    num//den ≡ a (mod m)

Writes the answer to `num` and `den` inplace. Returns `true` if the
reconstrction was successful and `false` otherwise. 

## Additional parameters:

- `bnd`: stores the stopping criterion threshold (see
    `rational_reconstruction_bound`) 
- `buf`, `buf1`, `buf2`, `buf3`, `u1`, `u2`, `u3`,  `v1`, `v2`, `v3`: buffers
"""
function rational_reconstruction!(
    num::BigInt,
    den::BigInt,
    bnd::BigInt,
    buf::BigInt,
    buf1::BigInt,
    buf2::BigInt,
    buf3::BigInt,
    u1::BigInt,
    u2::BigInt,
    u3::BigInt,
    v1::BigInt,
    v2::BigInt,
    v3::BigInt,
    a::BigInt,
    m::BigInt
)
    # Assumes the input is nonnegative!
    @invariant Base.GMP.MPZ.cmp_ui(a, 0) >= 0

    # fast path for numbers smaller than O(sqrt(modulo))
    if Base.GMP.MPZ.cmp(a, bnd) < 0
        Base.GMP.MPZ.set!(num, a)
        Base.GMP.MPZ.set_ui!(den, 1)
        return true
    end

    Base.GMP.MPZ.set_ui!(u1, 1)
    Base.GMP.MPZ.set_ui!(u2, 0)
    Base.GMP.MPZ.set!(u3, m)
    Base.GMP.MPZ.set_ui!(v1, 0)
    Base.GMP.MPZ.set_ui!(v2, 1)
    Base.GMP.MPZ.set!(v3, a)

    while true
        if Base.GMP.MPZ.cmp(v2, bnd) > 0
            return false
        end

        Base.GMP.MPZ.set!(buf, v3)
        if Base.GMP.MPZ.cmp_ui(buf, 0) < 0
            Base.GMP.MPZ.neg!(buf)
        end

        if Base.GMP.MPZ.cmp(buf, bnd) < 0
            break
        end

        Base.GMP.MPZ.tdiv_q!(buf, u3, v3)

        Base.GMP.MPZ.mul!(buf1, buf, v1)
        Base.GMP.MPZ.mul!(buf2, buf, v2)
        Base.GMP.MPZ.mul!(buf3, buf, v3)

        Base.GMP.MPZ.sub!(buf1, u1, buf1)
        Base.GMP.MPZ.sub!(buf2, u2, buf2)
        Base.GMP.MPZ.sub!(buf3, u3, buf3)

        Base.GMP.MPZ.set!(u1, v1)
        Base.GMP.MPZ.set!(u2, v2)
        Base.GMP.MPZ.set!(u3, v3)

        Base.GMP.MPZ.set!(v1, buf1)
        Base.GMP.MPZ.set!(v2, buf2)
        Base.GMP.MPZ.set!(v3, buf3)
    end

    Base.GMP.MPZ.set!(den, v2)
    Base.GMP.MPZ.set!(num, v3)

    #=
    Base.GMP.MPZ.gcd!(buf, den, num)
    Base.GMP.MPZ.tdiv_q!(den, buf)
    Base.GMP.MPZ.tdiv_q!(num, buf)
    =#

    if Base.GMP.MPZ.cmp_ui(den, 0) < 0
        Base.GMP.MPZ.neg!(den)
        Base.GMP.MPZ.neg!(num)
    end

    true
end

###
# CRT

"""
Implements the linear Chinese remainder lifting algorithm. Computes the unique
`x` such that
        
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
end

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
end

function CRT_precompute!(
    M::BigInt,
    n1::BigInt,
    n2::BigInt,
    ci::Vector{BigInt},
    moduli::Vector{UInt}
)
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
end
