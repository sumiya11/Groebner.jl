# This file is a part of Groebner.jl. License is GNU GPL v2.

###
# Backend for `normalform`

function normalform0(polynomials, to_be_reduced, options)
    ring, monoms, coeffs, options = io_convert_polynomials_to_ir(polynomials, options)
    # TODO: this deepcopy is sad
    ring_to_be_reduced, monoms_to_be_reduced, coeffs_to_be_reduced, _ =
        io_convert_polynomials_to_ir(to_be_reduced, deepcopy(options))
    monoms_reduced, coeffs_reduced = _normalform1(
        ring,
        monoms,
        coeffs,
        ring_to_be_reduced,
        monoms_to_be_reduced,
        coeffs_to_be_reduced,
        options
    )
    result = io_convert_ir_to_polynomials(
        ring,
        polynomials,
        monoms_reduced,
        coeffs_reduced,
        options
    )
    result
end

# (exponent vectors, coefficients) => (exponent vectors, coefficients)
function normalform1(
    ring,
    monoms,
    coeffs,
    ring_to_be_reduced,
    monoms_to_be_reduced,
    coeffs_to_be_reduced,
    options
)
    ring, monoms, coeffs = ir_ensure_assumptions(ring, monoms, coeffs)
    ring_to_be_reduced, monoms_to_be_reduced, coeffs_to_be_reduced = ir_ensure_assumptions(
        ring_to_be_reduced,
        monoms_to_be_reduced,
        coeffs_to_be_reduced
    )
    _normalform1(
        ring,
        monoms,
        coeffs,
        ring_to_be_reduced,
        monoms_to_be_reduced,
        coeffs_to_be_reduced,
        options
    )
end

function _normalform1(
    ring,
    monoms,
    coeffs,
    ring_to_be_reduced,
    monoms_to_be_reduced,
    coeffs_to_be_reduced,
    options
)
    try
        params = AlgorithmParameters(ring, options)
        return __normalform1(
            ring,
            monoms,
            coeffs,
            ring_to_be_reduced,
            monoms_to_be_reduced,
            coeffs_to_be_reduced,
            params
        )
    catch err
        if isa(err, MonomialDegreeOverflow)
            @log :info """
            Possible overflow of exponent vector detected. 
            Restarting with at least 32 bits per exponent."""
            params = AlgorithmParameters(ring, options; hint=:large_exponents)
            return __normalform1(
                ring,
                monoms,
                coeffs,
                ring_to_be_reduced,
                monoms_to_be_reduced,
                coeffs_to_be_reduced,
                params
            )
        else
            # Something bad happened.
            rethrow(err)
        end
    end
end

function __normalform1(
    ring,
    monoms,
    coeffs,
    ring_to_be_reduced,
    monoms_to_be_reduced,
    coeffs_to_be_reduced,
    params
)
    @invariant ir_is_valid(ring, monoms, coeffs)
    @invariant ir_is_valid(ring_to_be_reduced, monoms_to_be_reduced, coeffs_to_be_reduced)

    _, ring2, monoms2, coeffs2 =
        io_convert_ir_to_internal(ring, monoms, coeffs, params, params.representation)
    _, ring_to_be_reduced2, monoms_to_be_reduced2, coeffs_to_be_reduced2 =
        io_convert_ir_to_internal(
            ring_to_be_reduced,
            monoms_to_be_reduced,
            coeffs_to_be_reduced,
            params,
            params.representation
        )
    monoms_reduced2, coeffs_reduced2 = normalform2(
        ring2,
        monoms2,
        coeffs2,
        ring_to_be_reduced2,
        monoms_to_be_reduced2,
        coeffs_to_be_reduced2,
        params
    )
    monoms_reduced, coeffs_reduced = io_convert_internal_to_ir(
        ring_to_be_reduced2,
        monoms_reduced2,
        coeffs_reduced2,
        params
    )
    monoms_reduced, coeffs_reduced
end

function normalform2(
    ring,
    monoms,
    coeffs,
    ring_to_be_reduced,
    monoms_to_be_reduced,
    coeffs_to_be_reduced,
    params
)
    @invariant ring.nvars == ring_to_be_reduced.nvars &&
               ring.ch == ring_to_be_reduced.ch &&
               isequal(ring.ord, ring_to_be_reduced.ord)

    if params.check
        if !isgroebner2(ring, monoms, coeffs, params)
            throw(DomainError("Input polynomials do not look like a Groebner basis."))
        end
    end

    _monoms = filter(!isempty, monoms)
    _coeffs = filter(!isempty, coeffs)
    if isempty(_monoms)
        return monoms_to_be_reduced, coeffs_to_be_reduced
    end
    monoms, coeffs = _monoms, _coeffs

    nonzero_indices = findall(!io_iszero_coeffs, coeffs_to_be_reduced)
    if isempty(nonzero_indices)
        return monoms_to_be_reduced, coeffs_to_be_reduced
    end

    monoms_to_be_reduced_nonzero = monoms_to_be_reduced[nonzero_indices]
    coeffs_to_be_reduced_nonzero = coeffs_to_be_reduced[nonzero_indices]

    monoms_reduced, coeffs_reduced = _normalform2(
        ring,
        monoms,
        coeffs,
        monoms_to_be_reduced_nonzero,
        coeffs_to_be_reduced_nonzero,
        params
    )

    monoms_to_be_reduced[nonzero_indices] .= monoms_reduced
    coeffs_to_be_reduced[nonzero_indices] .= coeffs_reduced
    monoms_reduced = monoms_to_be_reduced
    coeffs_reduced = coeffs_to_be_reduced

    monoms_reduced, coeffs_reduced
end

function _normalform2(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    monoms_to_be_reduced::Vector{Vector{M}},
    coeffs_to_be_reduced::Vector{Vector{C}},
    params
) where {M <: Monom, C <: Coeff}
    basis, _, hashtable = f4_initialize_structs(ring, monoms, coeffs, params)
    tobereduced = basis_initialize_using_existing_hashtable(
        ring,
        monoms_to_be_reduced,
        coeffs_to_be_reduced,
        hashtable
    )
    f4_normalform!(ring, basis, tobereduced, hashtable, params.arithmetic)
    basis_export_data(tobereduced, hashtable)
end
