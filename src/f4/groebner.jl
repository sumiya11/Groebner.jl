
# Mutate everything!

#######################################
# Finite field groebner

function groebner_ff(
            ring::PolyRing,
            exps::Vector{Vector{Vector{UInt16}}},
            coeffs::Vector{Vector{UInt64}},
            reduced::Bool,
            rng::Rng) where {Rng<:Random.AbstractRNG}
    # specialize on ordering (not yet)
    # groebner_ff(ring, exps, coeffs, reduced, rng, Val(ring.ord))
    f4(ring, exps, coeffs, rng, reduced)
end

#######################################
# Rational field groebner

function groebner_qq(
            ring::PolyRing,
            exps::Vector{Vector{Vector{UInt16}}},
            coeffs::Vector{Vector{Rational{BigInt}}},
            reduced::Bool,
            randomized::Bool,
            rng::Rng,
            ) where {Rng<:Random.AbstractRNG}

    # we can mutate coeffs and exps here

    gbcoeffs_accum = Vector{Vector{BigInt}}(undef, 0)
    gbcoeffs_qq = Vector{Vector{Rational{BigInt}}}(undef, 0)

    # scale coefficients to integer ring inplace
    coeffs_zz = scale_denominators!(coeffs)

    modulo = BigInt(1)
    prime = nextluckyprime(coeffs_zz)

    ring_ff, coeffs_ff = reduce_modulo(coeffs_zz, ring, prime)

    # TODO: 2^16
    gens_ff, ht = initialize_structures(ring_ff, exps, coeffs_ff, rng, 2^16)

    i = 1
    while true
        @info "$i: selected lucky prime $prime"

        gb_ff, ht = f4!(ring_ff, gens_ff, ht, reduced)

        @info "CRT modulo ($modulo, $(ring_ff.ch))"
        reconstruct_crt!(gbcoeffs_accum, modulo, gb_ff.coeffs, ring_ff)

        @info "Reconstructing modulo $modulo"
        reconstruct_modulo!(gbcoeffs_qq, gbcoeffs_accum, modulo)

        if true # correctness_check()
            break
        end

        prime = nextluckyprime(coeffs_zz, prime)
        reduce_modulo!(coeffs_zz, ring, prime, ring_ff, coeffs_ff)
        reinitialize_structures!(gens_ff, ht, coeffs_ff)
    end

    gbexps, gbcoeffs_qq
end

#######################################
# Finite field isgroebner

# UWU!
function isgroebner_ff(
            ring::PolyRing,
            exps::Vector{Vector{Vector{UInt16}}},
            coeffs::Vector{Vector{UInt64}},
            rng)

    isgroebner_f4(ring, exps, coeffs, rng)
end

#######################################
# Rational field groebner

# TODO
function isgroebner_qq(
            ring::PolyRing,
            exps::Vector{Vector{Vector{UInt16}}},
            coeffs::Vector{Vector{Rational{BigInt}}},
            randomized::Bool,
            rng)

    if randomized
        coeffs_zz = scale_denominators!(coeffs)
        goodprime = nextgoodprime(coeffs_zz, Int[], 2^30 + 3)
        ring_ff, coeffs_ff = reduce_modulo(coeffs_zz, ring, goodprime)
        isgroebner_f4(ring_ff, exps, coeffs_ff, rng)
    else
        error("Sorry, not randomized version is not implemented yet.")
    end
end
