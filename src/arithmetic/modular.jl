
#=
    The file contains functions to operate modular reconstrction and CRT
=#

#------------------------------------------------------------------------------

const FIRST_COMPUTE_PRIME = 2^31 - 1
const FIRST_CHECK_PRIME   = 2^30 + 3

#------------------------------------------------------------------------------

# Rational number reconstruction implementation borrowed from CLUE
# and modified a bit to suit the 'Modern Computer Algebra' definitions
# Returns a rational r // h of QQ field in a canonical form such that
#   r // h ≡ a (mod m)
#
# let n = max( λ(a), λ(m) ) , where λ(x) is a number of bits for x
# O(n^2)
function rational_reconstruction(a::I, m::I) where {I<:Union{Int, BigInt}}
    a = mod(a, m)
    if a == 0 || m == 0
        return QQ(0, 1)
    end
    if m < 0
        m = -m
    end
    if a < 0
        a = m - a
    end
    if a == 1
        return QQ(1, 1)
    end
    bnd = sqrt(float(m) / 2)

    U = (I(1), I(0), m)
    V = (I(0), I(1), a)
    while abs(V[3]) >= bnd
        q = div(U[3], V[3])
        # TODO: use MutableArithmetics
        T = U .- q .* V
        U = V
        V = T
    end

    t = abs(V[2])
    r = V[3] * sign(V[2])
    # changed from `<= bnd` to `<= m / bnd`
    # we can speed up this !
    if t <= m / bnd && gcd(r, t) == 1
        return QQ(r, t)
    end

    # TODO: not needed
    throw(DomainError(
        :($a//$m), "rational reconstruction of $a (mod $m) does not exist"
    ))

    return QQ(0, 1)
end

rational_bound(x::BigInt) = ceil(BigInt, sqrt(float(x) / 2))

# 0 allocations!
#=
    Computes the rational reconstruction of `a` mod `m`.
    Namely, a pair of numbers num, den , such that
        num//den ≡ a (mod m)

    Additional params:
        bnd  = stores stopping criterion threshold
        buf, buf1, buf2, buf3  = buffers
        u1, u2, u3 = buffers
        v1, v2, v3 = buffers
=#
function rational_reconstruction!(
            num::BigInt, den::BigInt, bnd::BigInt, buf::BigInt,
            buf1::BigInt, buf2::BigInt, buf3::BigInt,
            u1::BigInt, u2::BigInt, u3::BigInt,
            v1::BigInt, v2::BigInt, v3::BigInt,
            a::BigInt, m::BigInt)

    if Base.GMP.MPZ.cmp_ui(a, 0) == 0
        Base.GMP.MPZ.set_ui!(num, 0)
        Base.GMP.MPZ.set_ui!(den, 1)
        return nothing
    end
    # TODO
    # assumes input is nonnegative
    @assert Base.GMP.MPZ.cmp_ui(a, 0) > 0

    if Base.GMP.MPZ.cmp_ui(a, 1) == 0
        Base.GMP.MPZ.set_ui!(num, 1)
        Base.GMP.MPZ.set_ui!(den, 1)
        return nothing
    end

    Base.GMP.MPZ.set_ui!(u1, 1)
    Base.GMP.MPZ.set_ui!(u2, 0)
    Base.GMP.MPZ.set!(u3, m)
    Base.GMP.MPZ.set_ui!(v1, 0)
    Base.GMP.MPZ.set_ui!(v2, 1)
    Base.GMP.MPZ.set!(v3, a)

    while true
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

    # TODO

    Base.GMP.MPZ.gcd!(buf, den, num)
    Base.GMP.MPZ.tdiv_q!(den, buf)
    Base.GMP.MPZ.tdiv_q!(num, buf)

    if Base.GMP.MPZ.cmp_ui(den, 0) < 0
        Base.GMP.MPZ.neg!(den)
        Base.GMP.MPZ.neg!(num)
    end

    nothing
end

#------------------------------------------------------------------------------

# 0 allocations!
#=
    Computes the unique x s.t
        x ≡ a1 mod m1
        x ≡ a2 mod m2

    Additional params:
        M   = m1*m2
        buf = a buffer
        n1  = a buffer
        n2  = a buffer
=#
function CRT!(
            M::BigInt, buf::BigInt, n1::BigInt, n2::BigInt,
            a1::BigInt, m1::BigInt, a2::UInt, m2::BigInt)

    Base.GMP.MPZ.gcdext!(buf, n1, n2, m1, m2)

    Base.GMP.MPZ.mul!(buf, m1, n1)
    Base.GMP.MPZ.mul_ui!(n1, buf, a2)

    Base.GMP.MPZ.mul!(buf, m2, n2)
    Base.GMP.MPZ.mul!(n2, buf, a1)

    Base.GMP.MPZ.add!(buf, n1, n2)
    Base.GMP.MPZ.fdiv_r!(buf, M)
end

#------------------------------------------------------------------------------