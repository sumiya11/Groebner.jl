# This file is a part of Groebner.jl. License is GNU GPL v2.

###
# Backend for auxiliary functions: lead

# polynomials => polynomials
function lead0(polynomials, options::KeywordArguments)
    ring, monoms, coeffs, options = io_convert_polynomials_to_ir(polynomials, options)
    lead_monoms, lead_coeffs = _lead1(ring, monoms, coeffs, options)
    result = io_convert_ir_to_polynomials(ring, polynomials, lead_monoms, lead_coeffs, options)
    result
end

# (exponent vectors, coefficients) => (exponent vectors, coefficients)
function lead1(
    ring::PolyRing,
    monoms::Vector{Vector{Vector{I}}},
    coeffs::Vector{Vector{C}},
    options::KeywordArguments
) where {I <: Integer, C <: Coeff}
    ring, monoms, coeffs = ir_ensure_valid(ring, monoms, coeffs)
    _lead1(ring, monoms, coeffs, options)
end

function _lead1(
    ring::PolyRing,
    monoms::Vector{Vector{Vector{I}}},
    coeffs::Vector{Vector{C}},
    options::KeywordArguments
) where {I <: Integer, C <: Coeff}
    @invariant ir_is_valid(ring, monoms, coeffs)
    params = AlgorithmParameters(ring, options; hint=:large_exponents)
    _, ring2, monoms2, coeffs2 = ir_convert_ir_to_internal(ring, monoms, coeffs, params)
    lead_monoms2, lead_coeffs2 = lead2(ring2, monoms2, coeffs2, params)
    lead_monoms, lead_coeffs = ir_convert_internal_to_ir(ring2, lead_monoms2, lead_coeffs2, params)
    lead_monoms, lead_coeffs
end

# internal structs => internal structs
function lead2(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: Coeff}
    lead_monoms = map(f -> isempty(f) ? f : [first(f)], monoms)
    lead_coeffs = map(f -> isempty(f) ? f : [first(f)], coeffs)
    lead_monoms, lead_coeffs
end
