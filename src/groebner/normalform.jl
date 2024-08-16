# This file is a part of Groebner.jl. License is GNU GPL v2.

###
# Backend for `normalform`

function _normalform0(polynomials, to_be_reduced, kws::KeywordArguments)
    polynomial_repr =
        io_select_polynomial_representation(polynomials, kws, hint=:large_exponents)
    ring, var_to_index1, monoms, coeffs =
        io_convert_to_internal(polynomial_repr, polynomials, kws)
    if isempty(monoms)
        @log :misc "Input basis consisting of zero polynomials only."
        return to_be_reduced
    end
    ring_to_be_reduced, var_to_index2, monoms_to_be_reduced, coeffs_to_be_reduced =
        io_convert_to_internal(polynomial_repr, to_be_reduced, kws, dropzeros=false)
    var_to_index = merge(var_to_index1, var_to_index2)
    nonzero_indices = findall(!io_iszero_coeffs, coeffs_to_be_reduced)
    if isempty(nonzero_indices)
        @log :misc "Polynomials to be reduced are all zero."
        return to_be_reduced
    end
    monoms_to_be_reduced_nonzero = monoms_to_be_reduced[nonzero_indices]
    coeffs_to_be_reduced_nonzero = coeffs_to_be_reduced[nonzero_indices]
    params = AlgorithmParameters(ring, polynomial_repr, kws)
    ring, _ = io_set_monomial_ordering!(ring, var_to_index, monoms, coeffs, params)
    ring_, _ = io_set_monomial_ordering!(
        ring_to_be_reduced,
        var_to_index,
        monoms_to_be_reduced_nonzero,
        coeffs_to_be_reduced_nonzero,
        params
    )
    @invariant ring.nvars == ring_to_be_reduced.nvars && ring.ch == ring_to_be_reduced.ch
    if kws.check
        @log :misc "As `check=true` was provided, checking that the given input is indeed a Groebner basis"
        if !_isgroebner1(ring, monoms, coeffs, params)
            __not_a_basis_error(
                polynomials,
                "Input polynomials do not look like a Groebner basis."
            )
        end
    end
    @log :misc """
      Finalized polynomial rings:
      Basis: $ring
      To be reduced: $ring_"""
    @invariant ring.nvars == ring_.nvars &&
               ring.ch == ring_.ch &&
               isequal(ring.ord, ring_.ord)
    monoms_reduced, coeffs_reduced = _normalform1(
        ring,
        monoms,
        coeffs,
        monoms_to_be_reduced_nonzero,
        coeffs_to_be_reduced_nonzero,
        params
    )
    @log :all """After normal form: 
    Monoms: $monoms_reduced 
    Coeffs: $coeffs_reduced"""
    monoms_to_be_reduced[nonzero_indices] .= monoms_reduced
    coeffs_to_be_reduced[nonzero_indices] .= coeffs_reduced
    monoms_reduced = monoms_to_be_reduced
    coeffs_reduced = coeffs_to_be_reduced
    # TODO: remove `to_be_reduced` from arguments here
    res = io_convert_to_output(ring, to_be_reduced, monoms_reduced, coeffs_reduced, params)
    res
end

@timeit function _normalform1(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    monoms_to_be_reduced::Vector{Vector{M}},
    coeffs_to_be_reduced::Vector{Vector{C}},
    params
) where {M <: Monom, C <: Coeff}
    @log :debug "Initializing structs for F4"
    basis, _, hashtable = f4_initialize_structs(ring, monoms, coeffs, params)
    tobereduced = basis_initialize_using_existing_hashtable(
                ring,
        monoms_to_be_reduced,
        coeffs_to_be_reduced,
        hashtable
    )
    @log :all """
      Polynomials in internal representation before reduction:
      Basis: $basis
      To be reduced: $tobereduced
      """
    f4_normalform!(ring, basis, tobereduced, hashtable, params.arithmetic)
    basis_export_data(tobereduced, hashtable)
end
