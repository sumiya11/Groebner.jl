
#------------------------------------------------------------------------------

function common_denominator(coeffs::Vector{Rational{BigInt}})
    den = BigInt(1)
    for c in coeffs
        # use buffer?
        Base.GMP.MPZ.lcm!(den, denominator(c))
    end
    den
end

function scale_denominators(coeffs::Vector{Rational{BigInt}})
    intcoeffs = [BigInt(0) for _ in 1:length(coeffs)]
    # lcm of map(denominator, coeffs[i])
    # so that den × coeffs[i] is integer
    den = common_denominator(coeffs)
    buf = BigInt(1)
    for i in 1:length(coeffs)
        t = numerator(coeffs[i])
        Base.GMP.MPZ.tdiv_q!(buf, den, denominator(coeffs[i]))
        Base.GMP.MPZ.mul!(intcoeffs[i], t, buf)
    end
    intcoeffs
end

# scales denominators inplace, returning new vector of integer coefficients
function scale_denominators(coeffs::Vector{Vector{Rational{BigInt}}})
    intcoeffs = Vector{Vector{BigInt}}(undef, length(coeffs))
    for i in 1:length(coeffs)
        intcoeffs[i] = scale_denominators(coeffs[i])
    end
    intcoeffs
end

#------------------------------------------------------------------------------

# coerce each coeff in `coeffs` into residuals modulo `prime`
function reduce_modulo(
        coeffs::Vector{Vector{BigInt}},
        ring::PolyRing, prime::Int)

    ffcoeffs = Vector{Vector{UInt64}}(undef, length(coeffs))
    for i in 1:length(coeffs)
        ffcoeffs[i] = Vector{UInt64}(undef, length(coeffs[i]))
        reduce_modulo!(coeffs[i], prime, ffcoeffs[i])
    end
    # why abstract?
    # TODO
    ring_ff = PolyRing(ring.nvars, ring.explen, ring.ord, prime, ring.origring)
    ring_ff, ffcoeffs
end

function reduce_modulo!(coeffs, prime, coeffs_ff)
    @assert length(coeffs) == length(coeffs_ff)

    p   = BigInt(prime)
    buf = BigInt(0)
    c   = BigInt(0)
    for i in 1:length(coeffs)
        Base.GMP.MPZ.set!(c, coeffs[i])
        if Base.GMP.MPZ.cmp_ui(c, 0) < 0
            Base.GMP.MPZ.fdiv_q!(buf, c, p)
            Base.GMP.MPZ.neg!(buf)
            Base.GMP.MPZ.mul_ui!(buf, prime)
            Base.GMP.MPZ.add!(c, buf)
        end
        @assert c >= 0
        Base.GMP.MPZ.tdiv_r!(buf, c, p)
        coeffs_ff[i] = UInt64(buf)
    end
end

function reduce_modulo!(
            coeffs_zz, prime,
            ring_ff, basis::Basis)

    for i in 1:length(coeffs_zz)
        reduce_modulo!(coeffs_zz[i], prime, basis.coeffs[i])
    end

    basis.ch   = prime
    ring_ff.ch = prime
end

#------------------------------------------------------------------------------

function isluckyprime(
        coeffs::Vector{Vector{BigInt}},
        prime::Int)
    buf = BigInt(0)
    p   = BigInt(prime)
    for poly in coeffs
        for c in poly
            if Base.GMP.MPZ.tdiv_r!(buf, c, p) == 0
                 return false
            end
        end
    end
    return true
end

nextluckyprime(coeffs::Vector{Vector{BigInt}}) = nextluckyprime(coeffs, 1)

function nextluckyprime(
        coeffs::Vector{Vector{BigInt}},
        prevprime::Int)

    if prevprime == 1
        candidate = FIRST_COMPUTE_PRIME
    else
        candidate = nextprime(prevprime + 1)
    end
    @assert candidate isa Int

    if !isluckyprime(coeffs, candidate)
        candidate = nextluckyprime(coeffs, candidate)
    end

    candidate
end

#------------------------------------------------------------------------------

function nextgoodprime(
        coeffs::Vector{Vector{BigInt}},
        moduli::Vector{Int})

    nextgoodprime(coeffs, moduli, FIRST_CHECK_PRIME-1)
end

function nextgoodprime(
        coeffs::Vector{Vector{BigInt}},
        moduli::Vector{Int},
        prevprime::Int)

    p = nextluckyprime(coeffs, prevprime)

    if p >= FIRST_COMPUTE_PRIME
        error("Too large coefficient encountered in the basis. This error will be fixed soon!")
    end
    
    if p in moduli
        return nextgoodprime(coeffs, moduli, p)
    end
    return p
end

#------------------------------------------------------------------------------

function reconstruct_crt!(gbcoeffs_accum, coeffs_ff, ring_ff)
    resize!(gbcoeffs_accum, length(coeffs_ff))
    for i in 1:length(coeffs_ff)
        gbcoeffs_accum[i] = [BigInt(0) for _ in 1:length(coeffs_ff[i])]
        for j in 1:length(coeffs_ff[i])
            cf = coeffs_ff[i][j]
            Base.GMP.MPZ.set_ui!(gbcoeffs_accum[i][j], cf)
        end
    end
end

function reconstruct_crt!(gbcoeffs_accum, modulo, coeffs_ff, ring_ff)
    if modulo == 1
        reconstruct_crt!(gbcoeffs_accum, coeffs_ff, ring_ff)
        Base.GMP.MPZ.mul_ui!(modulo, ring_ff.ch)
        return nothing
    end

    ch = ring_ff.ch
    bigch = BigInt(ch)

    buf = BigInt(0)
    n1, n2 = BigInt(0), BigInt(0)
    M = BigInt(0)
    Base.GMP.MPZ.mul_ui!(M, modulo, ch)

    for i in 1:length(coeffs_ff)
        for j in 1:length(coeffs_ff[i])
            ca = gbcoeffs_accum[i][j]
            cf = coeffs_ff[i][j]
            CRT!(M, buf, n1, n2, ca, modulo, cf, bigch)
            Base.GMP.MPZ.set!(gbcoeffs_accum[i][j], buf)
            @assert gbcoeffs_accum[i][j] >= 0
        end
    end
    Base.GMP.MPZ.mul_ui!(modulo, ring_ff.ch)

    return nothing
end

#------------------------------------------------------------------------------

function reconstruct_modulo!(gbcoeffs_qq, gbcoeffs_accum, modulo)
    if isempty(gbcoeffs_qq)
        resize!(gbcoeffs_qq, length(gbcoeffs_accum))
        for i in 1:length(gbcoeffs_accum)
            gbcoeffs_qq[i] = [Rational{BigInt}(0)
                                for _ in 1:length(gbcoeffs_accum[i])]
        end
    end

    bnd = rational_bound(modulo)
    buf, buf1, buf2, buf3  = (BigInt(0) for _ in 1:4)
    u1, u2, u3, v1, v2, v3 = (BigInt(0) for _ in 1:6)

    for i in 1:length(gbcoeffs_qq)
        for j in 1:length(gbcoeffs_qq[i])
            cz = gbcoeffs_accum[i][j]
            cq = gbcoeffs_qq[i][j]
            num, den = numerator(cq), denominator(cq)
            rational_reconstruction!(num, den, bnd, buf,
                                        buf1, buf2, buf3,
                                        u1, u2, u3, v1, v2, v3,
                                        cz, modulo)
            # TODO
            @assert gcd(numerator(cq), denominator(cq)) == 1
        end
    end

    nothing
end
