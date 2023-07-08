
function _normalform(polynomials, to_be_reduced, kws::KeywordsHandler)
    polynomial_repr =
        select_polynomial_representation(polynomials, kws, hint=:large_exponents)
    ring, monoms, coeffs = convert_to_internal(polynomial_repr, polynomials, kws)
    if isempty(monoms)
        @log "Input basis consisting of zero polynomials."
        return to_be_reduced
    end
    ring_to_be_reduced, monoms_to_be_reduced, coeffs_to_be_reduced =
        convert_to_internal(polynomial_repr, to_be_reduced, kws)
    nonzeroindices = findall(!iszero_coeffs, coeffs_to_be_reduced)
    if isempty(nonzeroindices)
        @log "Polynomials to be reduced are all zero."
        return to_be_reduced
    end
    monoms_to_be_reduced_nonzero = monoms_to_be_reduced[nonzeroindices]
    coeffs_to_be_reduced_nonzero = coeffs_to_be_reduced[nonzeroindices]
    params = AlgorithmParameters(ring, kws)
    ring = change_ordering_if_needed!(ring, monoms, coeffs, params)
    ring = change_ordering_if_needed!(
        ring_to_be_reduced,
        monoms_to_be_reduced_nonzero,
        coeffs_to_be_reduced_nonzero,
        params
    )
    if kws.check
        @log "Checking that the given input is a Groebner basis"
        @assert _isgroebner(ring, monoms, coeffs, params) "The given set of polynomials is not a Groebner basis"
    end
    # TODO: check that ring and ring_to_be_reduced agree!
    monoms_reduced, coeffs_reduced = _normalform(
        ring,
        monoms,
        coeffs,
        monoms_to_be_reduced,
        coeffs_to_be_reduced,
        params
    )
    monoms_to_be_reduced[nonzeroindices] .= monoms_reduced
    coeffs_to_be_reduced[nonzeroindices] .= coeffs_reduced
    monoms_reduced = monoms_to_be_reduced
    coeffs_reduced = coeffs_to_be_reduced
    # TODO: remove `tobereduced` from here
    convert_to_output(ring, to_be_reduced, monoms_reduced, coeffs_reduced, kws)
end

function _normalform(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    monoms_to_be_reduced::Vector{Vector{M}},
    coeffs_to_be_reduced::Vector{Vector{C}},
    params
) where {M <: Monom, C <: Coeff}
    @log "Initializing structs for F4"
    basis, _, hashtable = initialize_structs(ring, monoms, coeffs, params)
    tobereduced = initialize_basis_using_existing_hashtable(
        ring,
        monoms_to_be_reduced,
        coeffs_to_be_reduced,
        hashtable
    )
    f4_normalform!(ring, basis, tobereduced, hashtable)
    export_basis_data(tobereduced, hashtable)
end
