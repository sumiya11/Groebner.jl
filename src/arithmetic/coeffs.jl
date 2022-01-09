
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
        intcoeffs[i] = den * numerator(coeffs[i])
    end
    intcoeffs
end

# scales denominators inplace, returning new vector of integer coefficients
function scale_denominators!(coeffs::Vector{Vector{Rational{BigInt}}})
    intcoeffs = Vector{Vector{BigInt}}(undef, length(coeffs))
    for i in 1:length(coeffs)
        intcoeffs[i] = remove_content!(coeffs[i])
    end
    intcoeffs
end

#------------------------------------------------------------------------------

# coerce each coeff in `coeffs` into residuals modulo `prime`
function coerce_coeffs(coeffs::Vector{BigInt}, prime::Int)
    ffcoeffs = Vector{UInt64}(undef, length(coeffs))
    for i in 1:length(coeffs)
        # TODO: faster
        # convert is guaranteed to be exact
        ffcoeffs[i] = UInt64((coeffs[i] + prime) % prime)
    end
    ffcoeffs
end

# coerce each coeff in `coeffs` into residuals modulo `prime`
function reduce_modulo(
        coeffs::Vector{Vector{BigInt}},
        ring::PolyRing, prime::Int)

    ffcoeffs = Vector{Vector{UInt64}}(undef, length(coeffs))
    for i in 1:length(coeffs)
        ffcoeffs[i] = coerce_coeffs(coeffs[i], prime)
    end
    ring_ff = PolyRing(ring.nvars, ring.explen, ring.ord, prime)
    ring_ff, ffcoeffs
end

#------------------------------------------------------------------------------

# checks if coefficients in coeffs are not large compared to modulo
# Having relatively large coeffs would mean that reconstruction probably failed
function reconstruction_check(
        coeffs::Vector{Vector{Rational{BigInt}}},
        modulo::BigInt
    )

    @inbounds for i in 1:length(coeffs)
        for j in 1:length(coeffs[i])
            c = coeffs[i][j]
            pval = ( numerator(c) * denominator(c) ) ^ 2
            if pval >= modulo
                return false
            end
        end
    end

    return true
end
