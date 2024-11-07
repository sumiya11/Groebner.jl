# This file is a part of Groebner.jl. License is GNU GPL v2.

###
# Backend for `normalform`

# polynomials => polynomials
function normalform0(
    polynomials::AbstractVector,
    to_be_reduced::AbstractVector,
    options::KeywordArguments
)
    # TODO: this deepcopy is sad
    _options = deepcopy(options)
    ring, monoms, coeffs, options = io_convert_polynomials_to_ir(polynomials, options)
    ring_tbr, monoms_tbr, coeffs_tbr, _ = io_convert_polynomials_to_ir(to_be_reduced, _options)
    monoms_reduced, coeffs_reduced =
        _normalform1(ring, monoms, coeffs, ring_tbr, monoms_tbr, coeffs_tbr, options)
    result =
        io_convert_ir_to_polynomials(ring, polynomials, monoms_reduced, coeffs_reduced, options)
    result
end

# (exponent vectors, coefficients) => (exponent vectors, coefficients)
function normalform1(
    ring::PolyRing,
    monoms::Vector{Vector{Vector{I1}}},
    coeffs::Vector{Vector{C1}},
    ring_tbr::PolyRing,
    monoms_tbr::Vector{Vector{Vector{I2}}},
    coeffs_tbr::Vector{Vector{C2}},
    options::KeywordArguments
) where {I1 <: Integer, C1 <: Coeff, I2 <: Integer, C2 <: Coeff}
    ring, monoms, coeffs = ir_ensure_valid(ring, monoms, coeffs)
    ring_tbr, monoms_tbr, coeffs_tbr = ir_ensure_valid(ring_tbr, monoms_tbr, coeffs_tbr)
    _normalform1(ring, monoms, coeffs, ring_tbr, monoms_tbr, coeffs_tbr, options)
end

function _normalform1(
    ring::PolyRing,
    monoms::Vector{Vector{Vector{I1}}},
    coeffs::Vector{Vector{C1}},
    ring_tbr::PolyRing,
    monoms_tbr::Vector{Vector{Vector{I2}}},
    coeffs_tbr::Vector{Vector{C2}},
    options::KeywordArguments
) where {I1 <: Integer, C1 <: Coeff, I2 <: Integer, C2 <: Coeff}
    try
        params = AlgorithmParameters(ring, options)
        return __normalform1(ring, monoms, coeffs, ring_tbr, monoms_tbr, coeffs_tbr, params)
    catch err
        if isa(err, MonomialDegreeOverflow)
            @info """
            Possible overflow of exponent vector detected. 
            Restarting with at least 32 bits per exponent.""" maxlog = 1
            params = AlgorithmParameters(ring, options; hint=:large_exponents)
            return __normalform1(ring, monoms, coeffs, ring_tbr, monoms_tbr, coeffs_tbr, params)
        else
            # Something bad happened.
            rethrow(err)
        end
    end
end

function __normalform1(
    ring::PolyRing,
    monoms::Vector{Vector{Vector{I1}}},
    coeffs::Vector{Vector{C1}},
    ring_tbr::PolyRing,
    monoms_tbr::Vector{Vector{Vector{I2}}},
    coeffs_tbr::Vector{Vector{C2}},
    params::AlgorithmParameters
) where {I1 <: Integer, C1 <: Coeff, I2 <: Integer, C2 <: Coeff}
    @invariant ir_is_valid(ring, monoms, coeffs)
    @invariant ir_is_valid(ring_tbr, monoms_tbr, coeffs_tbr)

    _, ring2, monoms2, coeffs2 = ir_convert_ir_to_internal(ring, monoms, coeffs, params)
    _, ring_tbr2, monoms_tbr2, coeffs_tbr2 =
        ir_convert_ir_to_internal(ring_tbr, monoms_tbr, coeffs_tbr, params)
    monoms_reduced2, coeffs_reduced2 =
        normalform2(ring2, monoms2, coeffs2, ring_tbr2, monoms_tbr2, coeffs_tbr2, params)
    monoms_reduced, coeffs_reduced =
        ir_convert_internal_to_ir(ring_tbr2, monoms_reduced2, coeffs_reduced2, params)
    monoms_reduced, coeffs_reduced
end

# internal structs => internal structs
function normalform2(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    ring_tbr::PolyRing,
    monoms_tbr::Vector{Vector{M}},
    coeffs_tbr::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: Coeff}
    @invariant ring.nvars == ring_tbr.nvars &&
               ring.ch == ring_tbr.ch &&
               isequal(ring.ord, ring_tbr.ord)

    if params.check
        if !isgroebner2(ring, monoms, coeffs, params)
            throw(DomainError("Input polynomials do not look like a Groebner basis."))
        end
    end

    _monoms = filter(!isempty, monoms)
    _coeffs = filter(!isempty, coeffs)
    if isempty(_monoms)
        return monoms_tbr, coeffs_tbr
    end
    monoms, coeffs = _monoms, _coeffs

    nonzero_indices = findall(!isempty, coeffs_tbr)
    if isempty(nonzero_indices)
        return monoms_tbr, coeffs_tbr
    end

    monoms_tbr_nonzero = monoms_tbr[nonzero_indices]
    coeffs_tbr_nonzero = coeffs_tbr[nonzero_indices]

    monoms_reduced, coeffs_reduced =
        _normalform2(ring, monoms, coeffs, monoms_tbr_nonzero, coeffs_tbr_nonzero, params)

    monoms_tbr[nonzero_indices] .= monoms_reduced
    coeffs_tbr[nonzero_indices] .= coeffs_reduced
    monoms_reduced = monoms_tbr
    coeffs_reduced = coeffs_tbr

    monoms_reduced, coeffs_reduced
end

function _normalform2(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    monoms_tbr::Vector{Vector{M}},
    coeffs_tbr::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: Coeff}
    basis, _, hashtable = f4_initialize_structs(ring, monoms, coeffs, params)
    tobereduced = basis_initialize_using_existing_hashtable(ring, monoms_tbr, coeffs_tbr, hashtable)
    f4_normalform!(ring, basis, tobereduced, hashtable, params.arithmetic)
    basis_export_data(tobereduced, hashtable)
end
