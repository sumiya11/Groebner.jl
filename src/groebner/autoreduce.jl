# This file is a part of Groebner.jl. License is GNU GPL v2.

### 
# Backend for `autoreduce`

function autoreduce0(polynomials, options)
    ring, monoms, coeffs, options = io_convert_polynomials_to_ir(polynomials, options)
    res_monoms, res_coeffs = autoreduce1(ring, monoms, coeffs, options)
    result = io_convert_ir_to_polynomials(ring, polynomials, res_monoms, res_coeffs, options)
    result
end

function autoreduce1(ring, monoms, coeffs, options)
    ring, monoms, coeffs = ir_ensure_valid(ring, monoms, coeffs)
    _autoreduce1(ring, monoms, coeffs, options)
end

function _autoreduce1(ring, monoms, coeffs, options)
    try
        params = AlgorithmParameters(ring, options)
        return __autoreduce1(ring, monoms, coeffs, params)
    catch err
        if isa(err, MonomialDegreeOverflow)
            @info """
            Possible overflow of exponent vector detected. 
            Restarting with at least 32 bits per exponent.""" maxlog = 1
            params = AlgorithmParameters(ring, options; hint=:large_exponents)
            return __autoreduce1(ring, monoms, coeffs, params)
        else
            # Something bad happened.
            rethrow(err)
        end
    end
end

function __autoreduce1(ring, monoms, coeffs, params)
    @invariant ir_is_valid(ring, monoms, coeffs)
    _, ring2, monoms2, coeffs2 = ir_convert_ir_to_internal(ring, monoms, coeffs, params)
    res_monoms2, res_coeffs2 = autoreduce2(ring2, monoms2, coeffs2, params)
    res_monoms, res_coeffs = ir_convert_internal_to_ir(ring2, res_monoms2, res_coeffs2, params)
    res_monoms, res_coeffs
end

function autoreduce2(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: Coeff}
    _monoms = filter(!isempty, monoms)
    _coeffs = filter(!isempty, coeffs)
    if isempty(_monoms)
        return [monoms[1]], [coeffs[1]]
    end
    monoms, coeffs = _monoms, _coeffs
    _autoreduce2(ring, monoms, coeffs, params)
end

function _autoreduce2(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: Coeff}
    basis, _, hashtable = f4_initialize_structs(ring, monoms, coeffs, params)
    basis_update!(basis, hashtable)
    matrix = matrix_initialize(ring, C)
    symbol_ht = hashtable_initialize_secondary(hashtable)
    f4_autoreduce!(ring, basis, matrix, hashtable, symbol_ht, params)
    basis_standardize!(ring, basis, hashtable, params.arithmetic, params.changematrix)
    basis_export_data(basis, hashtable)
end
