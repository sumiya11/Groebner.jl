
#------------------------------------------------------------------------------

function common_denominator(coeffs::Vector{Rational{BigInt}})
    den = BigInt(1)
    for c in coeffs
        # TODO: Modular arithmetic
        den = lcm(den, denominator(c))
    end
    den
end

# given an array of rational coeffs
# scales them by the minimal number such that the coeffs are all integer
function remove_content!(coeffs::Vector{Rational{BigInt}})
    intcoeffs = Vector{BigInt}(undef, length(coeffs))
    # lcm of map(denominator, coeffs[i])
    # so that den × coeffs[i] is integer
    den = common_denominator(coeffs)
    for i in 1:length(coeffs)
        # TODO: MutableArithmetics
        delta = divexact(den, denominator(coeffs[i]))
        intcoeffs[i] = numerator(coeffs[i]) * delta
    end
    intcoeffs
end

function scale_denominators!(coeffs::Vector{Rational{BigInt}})
    remove_content!(coeffs)
end

# scales denominators inplace, returning new vector of integer coefficients
function scale_denominators!(coeffs::Vector{Vector{Rational{BigInt}}})
    intcoeffs = Vector{Vector{BigInt}}(undef, length(coeffs))
    for i in 1:length(coeffs)
        intcoeffs[i] = scale_denominators!(coeffs[i])
    end
    intcoeffs
end

#------------------------------------------------------------------------------

# coerce each coeff in `coeffs` into residuals modulo `prime`
function reduce_modulo(coeffs::Vector{BigInt}, prime::Int)
    ffcoeffs = Vector{UInt64}(undef, length(coeffs))
    for i in 1:length(coeffs)
        # TODO: faster
        # convert is guaranteed to be exact
        c = coeffs[i]

        if c < 0
            #println(c, " ", prime*(div(c, prime) + 1) )
            c = c - prime*(div(c, prime) - 1)
        end
        #println(c)
        ffcoeffs[i] = UInt64(c % prime)
    end
    ffcoeffs
end

# coerce each coeff in `coeffs` into residuals modulo `prime`
function reduce_modulo(
        coeffs::Vector{Vector{BigInt}},
        ring::PolyRing, prime::Int)

    ffcoeffs = Vector{Vector{UInt64}}(undef, length(coeffs))
    for i in 1:length(coeffs)
        ffcoeffs[i] = reduce_modulo(coeffs[i], prime)
    end
    ring_ff = PolyRing(ring.nvars, ring.explen, ring.ord, prime, :abstract)
    ring_ff, ffcoeffs
end

#------------------------------------------------------------------------------

# Not used
# checks if coefficients in coeffs are not large compared to modulo
# Having relatively large coeffs would mean that reconstruction probably failed
function reconstruction_check(
        coeffs::Vector{Vector{Rational{BigInt}}},
        modulo::BigInt
    )

    @inbounds for i in 1:length(coeffs)
        for j in 1:length(coeffs[i])
            c = coeffs[i][j]
            # heuristic
            # TODO: a better one
            pval = ( numerator(c) + denominator(c) ) ^ 2
            if pval >= modulo
                return false
            end
        end
    end

    return true
end

#------------------------------------------------------------------------------

function isluckyprime(
        coeffs::Vector{Vector{BigInt}},
        prime::Int)
    for poly in coeffs
        for c in poly
            if c % prime == 0
                 return false
            end
        end
    end
    return true
end

function nextluckyprime(
        coeffs::Vector{Vector{BigInt}},
        prevprime::Int)

    if prevprime == 1
        candidate = 2^31 - 1
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
        moduli::Vector{Int},
        prevprime::Int)

    p = nextluckyprime(coeffs, prevprime)
    if p in moduli
        return nextgoodprime(coeffs, moduli, p)
    end
    return p
end
