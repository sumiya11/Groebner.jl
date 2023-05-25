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

Writes the answer to `num` and `den` inplace.

Returns `true` if the reconstrction was successful and `false` otherwise. 

# Additional parameters:
- `bnd`: stores the stopping criterion threshold 
    (see `rational_reconstruction_bound`) 
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
    if Base.GMP.MPZ.cmp_ui(a, 0) == 0
        Base.GMP.MPZ.set_ui!(num, 0)
        Base.GMP.MPZ.set_ui!(den, 1)
        return true
    end

    # assumes the input is nonnegative
    @assert Base.GMP.MPZ.cmp_ui(a, 0) > 0

    # fast path for 1//1 if the input is 1
    if Base.GMP.MPZ.cmp_ui(a, 1) == 0
        Base.GMP.MPZ.set_ui!(num, 1)
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

"""
Implements the Chinese remainder algorithm.
Computes the unique `x` such that
        
    x ≡ a1 mod m1
    x ≡ a2 mod m2

Writes the answer to `buf` inplace.

# Additional params:
- `M`: must be equal `m1 * m2`
- `buf`, `n1`, `n2`: buffers
- `minv1`: must be equal `m1^-1 mod M`
- `minv2`: must be equal `m2^-1 mod M`
"""
function CRT!(
    M::BigInt,
    buf::BigInt,
    n1::BigInt,
    n2::BigInt,
    a1::BigInt,
    minv1::BigInt,
    a2::UInt,
    minv2::BigInt,
    m1::BigInt,
    m2::BigInt
)
    Base.GMP.MPZ.mul!(buf, m1, minv1)
    Base.GMP.MPZ.mul_ui!(n1, buf, a2)

    Base.GMP.MPZ.mul!(buf, m2, minv2)
    Base.GMP.MPZ.mul!(n2, buf, a1)

    Base.GMP.MPZ.add!(buf, n1, n2)
    Base.GMP.MPZ.fdiv_r!(buf, M)
end
