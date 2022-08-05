
#=
    The file contains rational reconstruction and CRT reconstruction functions
=#

#------------------------------------------------------------------------------

# A reference implementation of rational reconstruction
#
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

    @debug "" bnd
    U = (I(1), I(0), m)
    V = (I(0), I(1), a)
    while abs(V[3]) >= bnd
        q = div(U[3], V[3])
        T = U .- q .* V
        U = V
        V = T
    end

    t = abs(V[2])
    r = V[3] * sign(V[2])

    @debug "" r t

    # changed from `<= bnd` to `<= m / bnd`
    # we can speed up this !
    # if t <= bnd && gcd(r, t) == 1
    #    return QQ(r, t)
    # end

    return QQ(r, t)

    # throw(DomainError(
    #    :($a//$m), "rational reconstruction of $a (mod $m) does not exist"
    # ))

    # return QQ(0, 1)
end

# returns the bound for rational reconstruction
# based on the current modulo size
#
# The bound for rational reconstrction:
# as soon as the numerator in rational reconstruction
# exceeds this bound, the gcd iteration is stopped
function rational_reconstruction_bound(modulo::BigInt)
    setprecision(2*Base.GMP.MPZ.sizeinbase(modulo, 2)) do
        ceil(BigInt, sqrt(BigFloat(modulo) / 2))
    end
end

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
        return true
    end

    # assumes input is nonnegative
    @assert Base.GMP.MPZ.cmp_ui(a, 0) > 0

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
            # @error "" v2 bnd
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
    num and den are obviously coprime,
    gcd is not nessesary here

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

#------------------------------------------------------------------------------

# 0 allocations!
#=
    Computes the unique x s.t
        x ≡ a1 mod m1
        x ≡ a2 mod m2

    Additional params:
        M   = m1*m2
        buf = buffer
        n1  = buffer
        n2  = buffer
        minv1 = m1^-1 mod M
        minv2 = m2^-1 mod M
=#
function CRT!(
            M::BigInt, buf::BigInt, n1::BigInt, n2::BigInt,
            a1::BigInt, minv1::BigInt, a2::UInt, minv2::BigInt,
            m1::BigInt, m2::BigInt)

    Base.GMP.MPZ.mul!(buf, m1, minv1)
    Base.GMP.MPZ.mul_ui!(n1, buf, a2)

    Base.GMP.MPZ.mul!(buf, m2, minv2)
    Base.GMP.MPZ.mul!(n2, buf, a1)

    Base.GMP.MPZ.add!(buf, n1, n2)
    Base.GMP.MPZ.fdiv_r!(buf, M)
end

#------------------------------------------------------------------------------
