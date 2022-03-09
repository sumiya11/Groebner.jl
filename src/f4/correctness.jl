

function correctness_check!(exps, coeffs, coeffs_zz, ring_ff, gb_ff, initial_ff,
                                ht, gbcoeffs_qq, gbcoeffs_accum,
                                modulo, goodprime, meta, tablesize, rng)

    # @warn "" heuristic_correctness_check(gbcoeffs_qq, modulo)
    # @warn "" randomized_correctness_check!(deepcopy(coeffs_zz), ring_ff, deepcopy(gb_ff), deepcopy(initial_ff),ht, deepcopy(gbcoeffs_qq), goodprime)

    # global HEURISTIC
    # global RANDOMIZED

    # HEURISTIC += heuristic_correctness_check(gbcoeffs_qq, modulo)
    # RANDOMIZED += randomized_correctness_check!(deepcopy(coeffs_zz), ring_ff, deepcopy(gb_ff), deepcopy(initial_ff), ht, deepcopy(gbcoeffs_qq), goodprime)

    # first we check coefficients only
    if meta.heuristiccheck
        if !heuristic_correctness_check(gbcoeffs_qq, modulo)
            @info "Heuristic check failed."
            return false
        end
        @info "Heuristic check passed!"
    end

    if meta.randomizedcheck
        if !randomized_correctness_check!(coeffs_zz, ring_ff, gb_ff, initial_ff, ht, gbcoeffs_qq, goodprime)
            @info "Randomized check failed."
            return false
        end
        @info "Randomized check passed!"
    end

    #=
    # if !flag1 || !flag2 || !flag3
    if !randomized_correctness_check!(coeffs_zz, ring_ff, gb_ff, initial_ff, ht, gbcoeffs_qq, goodprime)
        @info "Randomized check failed."
        return false
    end
    @info "Randomized check passed!"

    if !heuristic_correctness_check(gbcoeffs_qq, modulo)
    # if !heuristic_correctness_check(gbcoeffs_qq, modulo)
        @info "Heuristic check failed."
        return false
    end
    @info "Heuristic check passed!"
    =#

    t = time()
    if meta.guaranteedcheck
        if guaranteed_correctness_check(gb_ff.gens, gbcoeffs_qq, exps, coeffs,
                                            ht, ring_ff, tablesize, rng)

            @info "Guaranteed check passed!"
            return true
        end


        @info "Guaranteed check failed."
        return false
    end

    return true
end

@inline function threshold_in_heuristic(sznum, szden, szmod)
    1.10*(sznum + sznum) >= szmod
end

# ln(num) + ln(den) < 2 ln p
function heuristic_correctness_check(gbcoeffs_qq, modulo)
    #=
        We take advantage of the fact that Base.GMP.MPZ.sizeinbase(x, 2)
        marks the position of the leading bit in x.
        So we do not need to compute logarithms explicitly
    =#

    lnm = Base.GMP.MPZ.sizeinbase(modulo, 2)
    for i in 1:length(gbcoeffs_qq)
        for j in 1:length(gbcoeffs_qq[i])
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

function randomized_correctness_check!(coeffs_zz, ring_ff, gb_ff, initial_ff,
                                            ht, gbcoeffs_qq, goodprime)

    # @warn "reducing in correctness" goodprime
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

    # @warn "inside correctness qq" gbcoeffs_qq
    # @warn "inside correctness zz" gbcoeffs_zz

    reduce_modulo!(gbcoeffs_zz, goodprime, ring_ff, gb_ff)

    # @warn "inside correctness ff" gb_ff.coeffs gb_ff.ch

    @assert ring_ff.ch == goodprime == gb_ff.ch

    #=
    @error "!!!"
    println(gb_ff)
    println(ht.exponents[1:10])
    =#

    normalize_basis!(gb_ff)
    normalize_basis!(initial_ff)

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

    # buf_ff = copy_basis(gb_ff)
    #=
    println(map(bitstring, gb_ff.lead))
    println(map(bitstring, buf_ff.lead))

    @error "hashtable condition"
    println(ht.exponents[2:ht.load])
    println(map(x->bitstring(x.divmask), ht.hashdata[2:ht.load]))
    =#

    # check that initial ideal contains in the computed groebner basis modulo goodprime
    normal_form_f4!(ring_ff, gb_ff, ht, initial_ff)

    #=
    println(gb_ff)
    println("###")
    println(buf_ff)

    println(map(bitstring, gb_ff.lead))
    println(map(bitstring, buf_ff.lead))
    =#

    # @error "normal form"
    # println(initial_ff)
    for i in 1:initial_ff.ndone
        # meaning that it is not reduced
        if !isempty(initial_ff.coeffs[i])
            return false
        end
    end

    @assert ring_ff.ch == gb_ff.ch

    # check that the basis is groebner basis modulo goodprime
    if !isgroebner_f4!(ring_ff, gb_ff, ht)
        return false
    end

    @assert ring_ff.ch == gb_ff.ch == goodprime

    return true
end

#------------------------------------------------------------------------------

function guaranteed_correctness_check(gbexps, gbcoeffs_qq,
                                exps, coeffs, ht, ring_qq, tablesize, rng)

    gens_qq, _ = initialize_structures(ring_qq, exps, coeffs, ht)
    gb_qq, _   = initialize_structures(ring_qq, gbexps, gbcoeffs_qq, ht)

    normalize_basis!(gb_qq)
    normalize_basis!(gens_qq)

    # check that initial ideal contains in the computed groebner basis modulo goodprime
    normal_form_f4!(ring_qq, gb_qq, ht, gens_qq)

    # @error "normal form"
    # println(initial_ff)
    for i in 1:gens_qq.ndone
        # meaning that it is not reduced
        if !isempty(gens_qq.coeffs[i])
            return false
        end
    end

    # check that the basis is groebner basis modulo goodprime
    if !isgroebner_f4!(ring_qq, gb_qq, ht)
        return false
    end

    return true
end
