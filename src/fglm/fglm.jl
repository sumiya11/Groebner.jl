# This file is a part of Groebner.jl. License is GNU GPL v2.

function _fglm0(polynomials, ordering_from, ordering_to, kws)
    ordering_from == ordering_to && return polynomials

    representation = io_select_polynomial_representation(polynomials, kws)
    ring, var_to_index, monoms, coeffs =
        io_convert_to_internal(representation, polynomials, kws)
    params =
        AlgorithmParameters(ring, representation, kws, orderings=(ring.ord, ordering_from))
    if isempty(monoms)
        @log :misc "Input consisting of zero polynomials."
        throw(DomainError("Input consisting of zero polynomials to Groebner.fglm."))
        return io_convert_to_output(ring, polynomials, monoms, coeffs, params)
    end
    if kws.check
        @log :misc "Checking if input is a Grobner basis"
        # TODO this is, perhaps, broken!
        if !isgroebner(polynomials, ordering=ordering_from, certify=false)
            throw(DomainError("Input is not a Groebner basis."))
        end
    end

    ring, _ = io_set_monomial_ordering!(ring, var_to_index, monoms, coeffs, params)
    params = AlgorithmParameters(
        ring,
        representation,
        kws,
        orderings=(ordering_from, ordering_to)
    )

    _ordering_to = io_convert_to_internal_monomial_ordering(var_to_index, ordering_to)
    new_monoms, new_coeffs = _fglm1(ring, monoms, coeffs, _ordering_to, params)

    res = io_convert_to_output(ring, polynomials, new_monoms, new_coeffs, params)

    res
end

function _fglm1(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    ordering_to,
    params::AlgorithmParameters
) where {M <: Monom, C <: Coeff}
    basis, _, ht = f4_initialize_structs(ring, monoms, coeffs, params)
    new_basis, _, new_ht = fglm_main!(ring, basis, ht, ordering_to, params)
    basis_export_data(new_basis, new_ht)
end

###
###

function _fglm_residuals_in_batch(polynomials, ordering_from, ordering_to, kws)
    ordering_from == ordering_to && return polynomials

    representation =
        io_select_polynomial_representation(polynomials, kws, hint=:large_exponents)
    ring, var_to_index, monoms, coeffs =
        io_convert_to_internal(representation, polynomials, kws)
    params =
        AlgorithmParameters(ring, representation, kws, orderings=(ring.ord, ordering_from))
    if isempty(monoms)
        @log :misc "Input consisting of zero polynomials."
        throw(DomainError("Input consisting of zero polynomials to Groebner.fglm."))
        return io_convert_to_output(ring, polynomials, monoms, coeffs, params)
    end
    if kws.check
        @log :misc "Checking if input is a Grobner basis"
        # TODO this is, perhaps, broken!
        if !isgroebner(polynomials, ordering=ordering_from, certify=false)
            throw(DomainError("Input is not a Groebner basis."))
        end
    end

    ring, _ = io_set_monomial_ordering!(ring, var_to_index, monoms, coeffs, params)
    params = AlgorithmParameters(
        ring,
        representation,
        kws,
        orderings=(ordering_from, ordering_to)
    )

    _ordering_to = io_convert_to_internal_monomial_ordering(var_to_index, ordering_to)
    new_monoms, new_coeffs =
        _fglm_residuals_in_batch_1(ring, monoms, coeffs, _ordering_to, params)

    res = io_convert_to_output(ring, polynomials, new_monoms, new_coeffs, params)

    res
end

@timeit function _fglm_residuals_in_batch_1(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    ordering_to,
    params::AlgorithmParameters
) where {M <: Monom, C <: Coeff}
    basis, _, ht = f4_initialize_structs(ring, monoms, coeffs, params)
    new_basis, _, new_ht = _fglm_residuals_in_batch!(ring, basis, ht, ordering_to, params)
    basis_export_data(new_basis, new_ht)
end
