

function correctness_check!(coeffs_zz, ring_ff, gb_ff, initial_ff, ht,
                            gbcoeffs_qq, modulo, randomized, goodprime)
    # first we check coefficients only
    if !heuristic_correctness_check(gbcoeffs_qq, modulo)
        @info "Heuristic check failed."
        return false
    end
    @info "Heuristic check passed!"

    # @warn "correctnes for " coeffs_zz gb_ff.coeffs initial_ff.coeffs

    #
    if !randomized_correctness_check!(coeffs_zz, ring_ff, gb_ff, initial_ff,
                                        ht, gbcoeffs_qq, goodprime)
        @info "Randomized check failed."
        return false
    end
    @info "Randomized check passed!"

    if !randomized
        # TODO
        return true
        if guaranteed_correctness_check(ring, exps, coeffs_zz, gbexps, gbcoeffs_qq, rng)
            @info "Proved check passed!"
            return true
        end
        @info "Proved check failed."
        return false
    end
    return true
end

# ln(num) + ln(den) < 2 ln p
function heuristic_correctness_check(coeffs_qq, modulo)
    #=
        We take advantage of the fact that Base.GMP.MPZ.sizeinbase(x, 2)
        marks the position of the leading bit in x.
        So we do not need to compute logarithms explicitly
    =#

    lnm = 2*Base.GMP.MPZ.sizeinbase(modulo, 2)
    for i in 1:length(coeffs_qq)
        for j in 1:length(coeffs_qq[i])
            c = coeffs_qq[i][j]
            num, den = numerator(c), denominator(c)
            if Base.GMP.MPZ.sizeinbase(num, 2) + Base.GMP.MPZ.sizeinbase(den, 2) > lnm
                return false
            end
        end
    end
    return true
end

function randomized_correctness_check!(coeffs_zz, ring_ff, gb_ff, initial_ff,
                                            ht, gbcoeffs_qq, goodprime)

    @warn "reducing in correctness" goodprime
    reduce_modulo!(coeffs_zz, goodprime, ring_ff, initial_ff)

    @assert ring_ff.ch == goodprime == initial_ff.ch

    # TODO

    #=
    @error "!"
    println(initial_ff)
    println(ht.exponents[1:10])
    # initial_ff.coeffs[1][1] = UInt(88)
    =#

    gbcoeffs_zz = scale_denominators(gbcoeffs_qq)

    @warn "inside correctness qq" gbcoeffs_qq
    @warn "inside correctness zz" gbcoeffs_zz

    reduce_modulo!(gbcoeffs_zz, goodprime, ring_ff, gb_ff)

    @warn "inside correctness ff" gb_ff.coeffs gb_ff.ch

    @assert ring_ff.ch == goodprime == gb_ff.ch

    #=
    @error "!!!"
    println(gb_ff)
    println(ht.exponents[1:10])
    =#

    normalize_basis!(gb_ff)

    #=
    println(gb_ff)
    println(ht.exponents[1:10])

    println("#########################")
    println("zz basis coeffs")
    for v in gbcoeffs_zz
        println(v)
    end
    println("####")
    println("ff basis coeffs")
    for v in gb_ff.coeffs
        println(v)
    end
    println("####")
    println("ff initial coeffs")
    for v in initial_ff.coeffs[1:initial_ff.ntotal]
        println(v)
    end
    println("####")
    println("#########################")
    =#

    # check that initial ideal contains in the computed groebner basis modulo goodprime
    normal_form_f4!(ring_ff, gb_ff, ht, initial_ff)
    # @error "normal form"
    # println(initial_ff)
    for i in 1:initial_ff.ndone
        # meaning that it is not reduced
        if !isempty(initial_ff.coeffs[i])
            return false
        end
    end

    # check that the basis is groebner basis modulo goodprime
    if !isgroebner_f4!(ring_ff, gb_ff, ht)
        return false
    end

    return true
end

#------------------------------------------------------------------------------

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
            gbcoeffs_qq, prevgoodprime, rng)

    goodprime = nextgoodprime(coeffs_zz, moduli, prevgoodprime)
    # @warn "goodprime" coeffs_zz moduli goodprime gbcoeffs_qq
    ring_ff, coeffs_ff = reduce_modulo(coeffs_zz, ring, goodprime)

    gbcoeffs_zz = scale_denominators(gbcoeffs_qq)

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
