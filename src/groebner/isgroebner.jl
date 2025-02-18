# This file is a part of Groebner.jl. License is GNU GPL v2.

### 
# Backend for `isgroebner`

function isgroebner0(polynomials, options)
    ring, monoms, coeffs, options = io_convert_polynomials_to_ir(polynomials, options)
    isgroebner1(ring, monoms, coeffs, options)
end

function isgroebner1(ring, monoms, coeffs, options)
    ring, monoms, coeffs = ir_ensure_valid(ring, monoms, coeffs)
    _isgroebner1(ring, monoms, coeffs, options)
end

function _isgroebner1(ring, monoms, coeffs, options)
    try
        params = AlgorithmParameters(ring, options)
        return __isgroebner1(ring, monoms, coeffs, params)
    catch err
        if isa(err, MonomialDegreeOverflow)
            @info """
            Possible overflow of exponent vector detected. 
            Restarting with at least 32 bits per exponent.""" maxlog = 1
            params = AlgorithmParameters(ring, options; hint=:large_exponents)
            return __isgroebner1(ring, monoms, coeffs, params)
        else
            # Something bad happened.
            rethrow(err)
        end
    end
end

function __isgroebner1(ring, monoms, coeffs, params)
    @invariant ir_is_valid(ring, monoms, coeffs)
    _, ring2, monoms2, coeffs2 = ir_convert_ir_to_internal(ring, monoms, coeffs, params)
    isgroebner2(ring2, monoms2, coeffs2, params)
end

function isgroebner2(ring, monoms, coeffs, params)
    _monoms = filter(!isempty, monoms)
    _coeffs = filter(!isempty, coeffs)
    isempty(_monoms) && return true
    monoms, coeffs = _monoms, _coeffs
    _isgroebner2(ring, monoms, coeffs, params)
end

# Finite fields
function _isgroebner2(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: CoeffZp}
    basis, pairset, hashtable = f4_initialize_structs(ring, monoms, coeffs, params)
    res = f4_isgroebner!(ring, basis, pairset, hashtable, params.arithmetic)
    res
end

# Generic fields
function _isgroebner2(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: CoeffGeneric}
    basis, pairset, hashtable = f4_initialize_structs(ring, monoms, coeffs, params)
    res = f4_isgroebner!(ring, basis, pairset, hashtable, params.arithmetic)
    res
end

# Rational numbers
function _isgroebner2(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: CoeffQQ}
    basis, pairset, hashtable = f4_initialize_structs(ring, monoms, coeffs, params)
    # If an honest computation over the rationals is requested
    if params.certify_check
        flag = f4_isgroebner!(ring, basis, pairset, hashtable, params.arithmetic)
        return flag
    end
    # Otherwise, check modulo primes
    basis_zz = clear_denominators!(basis, deepcopy=false)
    state = ModularState{BigInt, Rational{BigInt}, UInt32}(basis_zz.coeffs)
    prime = modular_random_prime(state, params.rng)
    ring_ff, basis_ff = modular_reduce_mod_p!(ring, basis_zz, prime, deepcopy=true)
    arithmetic = select_arithmetic(CoeffModular, prime, :auto, false)
    flag1 = f4_isgroebner!(ring_ff, basis_ff, pairset, hashtable, arithmetic)
    prime = modular_random_prime(state, params.rng)
    ring_ff, basis_ff = modular_reduce_mod_p!(ring, basis_zz, prime, deepcopy=true)
    arithmetic = select_arithmetic(CoeffModular, prime, :auto, false)
    flag2 = f4_isgroebner!(ring_ff, basis_ff, pairset, hashtable, arithmetic)
    if flag1 == flag2
        return flag1
    end
    prime = modular_random_prime(state, params.rng)
    ring_ff, basis_ff = modular_reduce_mod_p!(ring, basis_zz, prime, deepcopy=true)
    arithmetic = select_arithmetic(CoeffModular, prime, :auto, false)
    flag3 = f4_isgroebner!(ring_ff, basis_ff, pairset, hashtable, arithmetic)
    if flag1 == flag3
        return flag1
    else
        return flag2
    end
    flag1
end
