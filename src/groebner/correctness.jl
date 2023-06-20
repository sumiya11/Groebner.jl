
# If one of the primes in the modular computation 
# was unlucky all along
function _monte_carlo_error(msg)
    @warn msg
    throw(RecoverableException(msg))
end

function _not_a_basis_error(basis, msg)
    throw(DomainError(basis, msg))
end

# Checks if basis is groebner basis with a randomized algorithm
function _check_isgroebner(basis)
    !isgroebner(basis, certify=false) &&
        _not_a_basis_error(basis, "Input does not look like a groebner basis.")
end

# Checks that all computed bases are of the same shape,
# if that is not the case, throws
function basis_shape_control(gb_coeffs1::AbstractVector, gb_coeffs2::AbstractVector)
    if length(gb_coeffs1) != length(gb_coeffs2)
        _monte_carlo_error("Unlucky reduction to zero in probabilistic linear algebra.")
    end
    @inbounds for p in 1:length(gb_coeffs1)
        if length(gb_coeffs1[p]) != length(gb_coeffs2[p])
            _monte_carlo_error("Unlucky cancellation of basis coefficients modulo a prime.")
        end
    end
    nothing
end

# Heuristic check for correctness of the reconstructed basis over rationals
function rational_correctness_bound(modulo::BigInt)
    setprecision(2 * Base.GMP.MPZ.sizeinbase(modulo, 2)) do
        ceil(BigInt, BigFloat(modulo)^(1 / 10))
    end
end

function basis_coeff_control(primetracker, coeffaccum)
    modulo = primetracker.modulo
    bound = rational_correctness_bound(modulo)
    for c in coeffaccum.gb_coeffs_qq
        for cc in c
            if cc > bound
                return nothing
            end
        end
    end
    _monte_carlo_error("Unlucky cancellation of basis coefficients modulo a prime.")
    nothing
end

# Check that the basis is reconstructed correctly.
# There are 3 levels of checks:
#   - heuristic check (discards obviously bad cases),
#   - randomized check,
#   - certification.
#
# By default, only the first two are active, which gives the correct basis
# with a high probability
function correctness_check!(
    coeffaccum,
    coeffbuffer,
    primetracker,
    meta,
    ring,
    exps,
    coeffs,
    coeffs_zz,
    gens_temp_ff,
    gb_ff,
    ht
)

    # first we check coefficients with a heuristic check only
    if meta.heuristiccheck
        if !heuristic_correctness_check(coeffaccum.gb_coeffs_qq, primetracker.modulo)
            @info "Heuristic check failed."
            return false
        end
        @info "Heuristic check passed!"
    end

    # then check that a basis is also a basis modulo a prime
    if meta.randomizedcheck
        if !randomized_correctness_check!(
            coeffbuffer,
            coeffaccum,
            ring,
            coeffs_zz,
            gens_temp_ff,
            gb_ff,
            primetracker,
            ht
        )
            @info "Randomized check failed."
            return false
        end
        @info "Randomized check passed!"
    end

    if meta.guaranteedcheck
        return guaranteed_correctness_check!(
            ring,
            gb_ff.monoms,
            coeffaccum.gb_coeffs_qq,
            exps,
            coeffs,
            gens_temp_ff,
            ht
        )
    end

    return true
end

threshold_in_heuristic(sznum, szden, szmod) = 1.30 * (sznum + sznum) >= szmod

# Checks that 
# ln(num) + ln(den) < c ln(modulo)
# for all coefficients num/den
function heuristic_correctness_check(gbcoeffs_qq, modulo)
    lnm = Base.GMP.MPZ.sizeinbase(modulo, 2)
    @inbounds for i in 1:length(gbcoeffs_qq)
        for j in 1:length(gbcoeffs_qq[i])
            n = numerator(gbcoeffs_qq[i][j])
            d = denominator(gbcoeffs_qq[i][j])
            if threshold_in_heuristic(
                Base.GMP.MPZ.sizeinbase(n, 2),
                Base.GMP.MPZ.sizeinbase(d, 2),
                lnm
            )
                return false
            end
        end
    end
    return true
end

function randomized_correctness_check!(
    coeffbuffer,
    coeffaccum,
    ring,
    coeffs_zz,
    gens_temp_ff,
    gb_ff,
    primetracker,
    ht
)
    goodprime = nextgoodprime!(primetracker)

    reduce_modulo!(coeffbuffer, coeffs_zz, gens_temp_ff.coeffs, goodprime)
    gens_ff_copy = deepcopy_basis(gens_temp_ff)
    cleanup_basis!(ring, gens_ff_copy, goodprime)

    gb_coeffs_zz = scale_denominators(coeffbuffer, coeffaccum.gb_coeffs_qq)
    gb_ff_copy = deepcopy_basis(gb_ff)
    reduce_modulo!(coeffbuffer, gb_coeffs_zz, gb_ff_copy.coeffs, goodprime)
    cleanup_basis!(ring, gb_ff_copy, goodprime)

    # check that initial ideal contains in the computed groebner basis modulo goodprime
    normal_form_f4!(ring, gb_ff_copy, ht, gens_ff_copy)
    for i in 1:(gens_ff_copy.ndone)
        # meaning that something is not reduced
        if !iszero_coeffvector(gens_ff_copy.coeffs[i])
            return false
        end
    end

    # check that the basis is a groebner basis modulo goodprime
    if !isgroebner_f4!(ring, gb_ff_copy, ht)
        return false
    end

    return true
end

function guaranteed_correctness_check!(
    ring,
    gbexps,
    gb_coeffs_qq,
    exps,
    coeffs,
    gens_tmp_ff,
    ht
)
    @info "Setting parameter certify=true in groebner is not recommended."

    gens_qq, _ = initialize_structures(ring, gens_tmp_ff.monoms[1:(gens_tmp_ff.ntotal)], coeffs, ht)
    gb_qq, _   = initialize_structures(ring, gbexps, gb_coeffs_qq, ht)

    normalize_basis!(ring, gb_qq)
    normalize_basis!(ring, gens_qq)

    gens_qq_copy = deepcopy_basis(gens_qq)

    normal_form_f4!(ring, gb_qq, ht, gens_qq_copy)
    for i in 1:(gens_qq_copy.ndone)
        # meaning that it is not reduced
        if !iszero_coeffvector(gens_qq_copy.coeffs[i])
            return false
        end
    end

    if !isgroebner_f4!(ring, gb_qq, ht)
        return false
    end

    return true
end
