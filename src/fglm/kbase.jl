# This file is a part of Groebner.jl. License is GNU GPL v2.

function _kbase0(polynomials, kws)
    representation =
        io_select_polynomial_representation(polynomials, kws, hint=:large_exponents)
    ring, var_to_index, monoms, coeffs =
        io_convert_to_internal(representation, polynomials, kws)
    params = AlgorithmParameters(ring, representation, kws)
    if isempty(monoms)
        @log level = -2 "Input consisting of zero polynomials."
        throw(DomainError("Input consisting of zero polynomials to Groebner.kbase."))
        return io_convert_to_output(ring, polynomials, monoms, coeffs, params)
    end
    if kws.check
        @log level = -2 "Checking if a Grobner basis"
        if !isgroebner(polynomials)
            throw(DomainError("Input is not a Groebner basis."))
        end
    end

    ring, _ = io_set_monomial_ordering!(ring, var_to_index, monoms, coeffs, params)
    m, c = _kbase1(ring, monoms, coeffs, params)

    res = io_convert_to_output(ring, polynomials, m, c, params)

    res
end

@timeit function _kbase1(
    ring::PolyRing,
    basisexps::Vector{Vector{M}},
    basiscoeffs::Vector{Vector{C}},
    params
) where {M, C <: Coeff}
    basis, _, ht = f4_initialize_structs(ring, basisexps, basiscoeffs, params)
    basis, linbasis, ht = fglm_main!(ring, basis, ht, ring.ord, params)
    basis_export_data(linbasis, ht)
end
