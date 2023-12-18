
function _kbase(polynomials, kws)
    representation =
        select_polynomial_representation(polynomials, kws, hint=:large_exponents)
    ring, var_to_index, monoms, coeffs =
        convert_to_internal(representation, polynomials, kws)
    params = AlgorithmParameters(ring, representation, kws)
    if isempty(monoms)
        @log level = -2 "Input consisting of zero polynomials."
        throw(DomainError("Input consisting of zero polynomials to Groebner.kbase."))
        return convert_to_output(ring, polynomials, monoms, coeffs, params)
    end
    if kws.check
        @log level = -2 "Checking if a Grobner basis"
        if !isgroebner(polynomials)
            throw(DomainError("Input is not a Groebner basis."))
        end
    end
    ring, _ = set_monomial_ordering!(ring, var_to_index, monoms, coeffs, params)
    m, c = kbase_f4(ring, monoms, coeffs, params)
    res = convert_to_output(ring, polynomials, m, c, params)
    print_performance_counters(params.statistics)
    res
end

@timeit function kbase_f4(
    ring::PolyRing,
    basisexps::Vector{Vector{M}},
    basiscoeffs::Vector{Vector{C}},
    params
) where {M, C <: Coeff}
    basis, pairset, ht = f4_initialize_structs(ring, basisexps, basiscoeffs, params)
    basis, linbasis, ht = fglm_f4!(ring, basis, ht, ring.ord, params)
    basis_export_data(linbasis, ht)
end
