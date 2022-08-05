

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

    if meta.randomizedcheck
        if !randomized_correctness_check!(coeffbuffer, coeffaccum, ring,
                                            coeffs_zz, gens_temp_ff, gb_ff, primetracker, ht)
            @info "Randomized check failed."
            return false
        end
        @info "Randomized check passed!"
    end

    if meta.guaranteedcheck
        return guaranteed_correctness_check!(ring, gb_ff.gens, coeffaccum.gb_coeffs_qq,
                                            exps, coeffs, gens_temp_ff, ht)
    end

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
            ring, coeffs_zz, gens_temp_ff, gb_ff, primetracker, ht)

    goodprime = nextgoodprime!(primetracker)

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

function crt_issame_check!(coeffbuffer, coeffaccum, gbcoeffs_ff,
                                        modulo, prime)
    prevcoeffs = coeffaccum.prev_gb_coeffs_zz
    currcoeffs = coeffaccum.gb_coeffs_zz
    prevmodulo = div(modulo, prime)
    println("In same check")
    println(prevcoeffs)
    println(prevmodulo)
    println(currcoeffs)
    println(modulo)
    for i in 1:length(prevcoeffs)
        @inbounds for j in 1:length(prevcoeffs[i])
            if currcoeffs[i][j] != prevcoeffs[i][j]
                if (modulo - currcoeffs[i][j]) != (prevmodulo - prevcoeffs[i][j])
                    # if (currcoeffs[i][j] % prevmodulo) != prevcoeffs[i][j]
                        return false
                    # end
                end
            end
        end
    end
    return true
end

#------------------------------------------------------------------------------

function guaranteed_correctness_check!(ring, gbexps, gb_coeffs_qq,
                                        exps, coeffs, gens_tmp_ff, ht)

    @warn "Setting parameter certify=true in groebner is not recommended."

    #=
    println("GB\n", gbexps, "\n", gb_coeffs_qq)
    println("GENS\n", gens_tmp_ff, "\n", coeffs)

    println("HT $(gens_tmp_ff.ntotal)")
    println(ht.exponents[1:3])
    =#

    # tmp_exps = [copy(init_gens.gens[i]) for i in 1:init_gens.ntotal]
    # tmp_coeffs = [
    #                copy(coeffs[i])
    #                for i in 1:init_gens.ntotal]

    gens_qq, _ = initialize_structures(ring, gens_tmp_ff.gens[1:gens_tmp_ff.ntotal], coeffs, ht)
    gb_qq, _   = initialize_structures(ring, gbexps, gb_coeffs_qq, ht)

    normalize_basis!(gb_qq)
    normalize_basis!(gens_qq)

    gens_qq_copy = copy_basis_thorough(gens_qq)

    # check that initial ideal contains in the computed groebner basis modulo goodprime
    normal_form_f4!(ring, gb_qq, ht, gens_qq_copy)

    # TODO

    # @error "normal form"
    # println(initial_ff)
    for i in 1:gens_qq_copy.ndone
        # meaning that it is not reduced
        if !isempty(gens_qq_copy.coeffs[i])
            return false
        end
    end

    # check that the basis is groebner basis modulo goodprime
    if !isgroebner_f4!(ring, gb_qq, ht)
        return false
    end

    return true
end

#------------------------------------------------------------------------------

function linear_algebra_error(reason)
    throw(AssertionError("Probabilistic part of the algorithm failed, sorry. Please, submit a github issue. You can try to rerun with another random seed.\nReason estimated: $reason"))
end

function shape_control(gb_coeffs1, gb_coeffs2)
    if length(gb_coeffs1) != length(gb_coeffs2)
        linear_algebra_error("Unlucky reduction to zero in probabilistic linear algebra.")
    end

    for p in 1:length(gb_coeffs1)
        if length(gb_coeffs1[p]) != length(gb_coeffs2[p])
            linear_algebra_error("Unlucky cancellation of basis coefficients.")
        end
    end

    nothing
end

function rational_correctness_bound(modulo::BigInt)
    setprecision(2*Base.GMP.MPZ.sizeinbase(modulo, 2)) do
        ceil(BigInt, BigFloat(modulo)^(1/10))
    end
end

function coeff_control(primetracker, coeffaccum)
    modulo = primetracker.modulo
    bound = rational_correctness_bound(modulo)
    for c in coeffaccum.gb_coeffs_qq
        for cc in c
            if cc > bound
                return true
            end
        end
    end

    linear_algebra_error("Unlucky reduction to zero in probabilistic linear algebra.")

    return false
end
