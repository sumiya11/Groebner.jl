

function correctness_check!(coeffaccum, coeffbuffer, primetracker, meta,
                                ring, exps, coeffs, coeffs_zz, gens_temp_ff, gb_ff, ht)

    # first we check coefficients only
    if meta.heuristiccheck
        if !heuristic_correctness_check(coeffaccum.gb_coeffs_qq, primetracker.modulo)
            @info "Heuristic check failed."
            return false
        end
        @info "Heuristic check passed!"
    end

    goodprime = nextgoodprime!(primetracker)
    if meta.randomizedcheck
        if !randomized_correctness_check!(coeffbuffer, coeffaccum, ring,
                                            coeffs_zz, gens_temp_ff, gb_ff, goodprime, ht)
            @info "Randomized check failed."
            return false
        end
        @info "Randomized check passed!"
    end

    # return guaranteed_correctness_check(gbexps, coeffaccum.gb_coeffs_qq,
    #                                exps, coeffs, ht)

    return true
end

@inline function threshold_in_heuristic(sznum, szden, szmod)
    # coefficient 1.10
    1.30*(sznum + sznum) >= szmod
end

# ln(num) + ln(den) < c ln p
function heuristic_correctness_check(gbcoeffs_qq, modulo)
    #=

    =#

    lnm = Base.GMP.MPZ.sizeinbase(modulo, 2)
    for i in 1:length(gbcoeffs_qq)
        @inbounds for j in 1:length(gbcoeffs_qq[i])
            n = numerator(gbcoeffs_qq[i][j])
            d = denominator(gbcoeffs_qq[i][j])
            # println(Base.GMP.MPZ.sizeinbase(c, 2), " ", lnm)
            if threshold_in_heuristic(Base.GMP.MPZ.sizeinbase(n, 2), Base.GMP.MPZ.sizeinbase(d, 2), lnm)
                return false
            end
        end
    end
    return true
end

# TODO: cleanup
function randomized_correctness_check!(
            coeffbuffer, coeffaccum,
            ring, coeffs_zz, gens_temp_ff, gb_ff, goodprime, ht)

    reduce_modulo!(coeffbuffer, coeffs_zz, gens_temp_ff.coeffs, goodprime)
    gens_ff_copy = copy_basis_thorough(gens_temp_ff)
    cleanup_gens!(ring, gens_ff_copy, goodprime)

    # TODO: can encapsulate BigInt coefficients into coeffaccum
    gb_coeffs_zz = scale_denominators(coeffbuffer, coeffaccum.gb_coeffs_qq)
    gb_ff_copy = copy_basis_thorough(gb_ff)
    reduce_modulo!(coeffbuffer, gb_coeffs_zz, gb_ff_copy.coeffs, goodprime)
    cleanup_gens!(ring, gb_ff_copy, goodprime)

    @assert ring.ch == gb_ff_copy.ch == gens_ff_copy.ch

    # check that initial ideal contains in the computed groebner basis modulo goodprime
    normal_form_f4!(ring, gb_ff_copy, ht, gens_ff_copy)

    for i in 1:gens_ff_copy.ndone
        # meaning that it is not reduced
        if !isempty(gens_ff_copy.coeffs[i])
            return false
        end
    end

    # check that the basis is groebner basis modulo goodprime
    if !isgroebner_f4!(ring, gb_ff_copy, ht)
        return false
    end

    return true
end

#------------------------------------------------------------------------------

function guaranteed_correctness_check(gbexps, gb_coeffs_qq,
                                        exps, coeffs, ht)


    tmp_exps = [copy(init_gens.gens[i]) for i in 1:init_gens.ntotal]
    # tmp_coeffs = [
    #                copy(coeffs[i])
    #                for i in 1:init_gens.ntotal]

    gens_qq, _ = initialize_structures(ring_qq, tmp_exps, tmp_coeffs, ht)
    gb_qq, _   = initialize_structures(ring_qq, gbexps, gb_coeffs_qq, ht)

    normalize_basis!(gb_qq)
    normalize_basis!(gens_qq)

    # check that initial ideal contains in the computed groebner basis modulo goodprime
    normal_form_f4!(ring_qq, gb_qq, ht, gens_qq)

    # TODO

    # @error "normal form"
    # println(initial_ff)
    for i in 1:gens_qq.ndone
        # meaning that it is not reduced
        if !isempty(gens_qq.coeffs[i])
            return false
        end
    end

    gb_qq, _   = initialize_structures(ring_qq, gbexps, gbcoeffs_qq, ht)

    # check that the basis is groebner basis modulo goodprime
    if !isgroebner_f4!(ring_qq, gb_qq, ht)
        return false
    end

    return true
end
