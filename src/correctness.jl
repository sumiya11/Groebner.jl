

function correctness_check(
            ring, exps, coeffs_zz, gbexps, gbcoeffs_qq,
            moduli, randomized, rng)
    if randomized_correctness_check(ring, exps, coeffs_zz, gbexps, gbcoeffs_qq, moduli, rng)
        @info "Randomized check passed!"
        # if guaranteed checks are required
        if !randomized
            if guaranteed_correctness_check(ring, exps, coeffs_zz, gbexps, gbcoeffs_qq, moduli, rng)
                @info "Guaranteed checks pass!"
                return true
            end
            @info "Guaranteed checks didn't pass."
            return false
        else
            return true
        end
    end
    @info "Randomized check failed."
    return false
end

function randomized_correctness_check(
            ring, exps, coeffs_zz, gbexps,
            gbcoeffs_qq, moduli, rng)

    goodprime = nextgoodprime(coeffs_zz, moduli, 2^30)
    # @warn "goodprime" coeffs_zz moduli goodprime gbcoeffs_qq
    ring_ff, coeffs_ff = reduce_modulo(coeffs_zz, ring, goodprime)

    gbcoeffs_zz = scale_denominators!(gbcoeffs_qq)

    # @warn "gbzz" gbcoeffs_zz

    ring_ff, gbcoeffs_ff = reduce_modulo(gbcoeffs_zz, ring, goodprime)

    # check that initial ideal contains in the computed groebner basis modulo goodprime
    expsreduced, _ = normal_form_f4(ring_ff, gbexps, gbcoeffs_ff, exps, coeffs_ff, rng)
    for es in expsreduced
        if !isempty(es)
            return false
        end
    end

    # check that the basis is groebner basis modulo goodprime
    if !isgroebner_f4(ring_ff, gbexps, gbcoeffs_ff, rng)
        return false
    end

    return true
end

function guaranteed_correctness_check(
            ring, exps, coeffs_zz, gbexps,
            gbcoeffs_qq, moduli, rng)

end
