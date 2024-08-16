# This file is a part of Groebner.jl. License is GNU GPL v2.

### 
# Backend for `isgroebner`

function _isgroebner0(polynomials, kws::KeywordArguments)
    polynomial_repr =
        io_select_polynomial_representation(polynomials, kws, hint=:large_exponents)
    ring, var_to_index, monoms, coeffs =
        io_convert_to_internal(polynomial_repr, polynomials, kws)
    if isempty(monoms)
        @log :misc "Input consisting of zero polynomials, which is a Groebner basis by our convention"
        return true
    end
    params = AlgorithmParameters(ring, polynomial_repr, kws)
    ring, _ = io_set_monomial_ordering!(ring, var_to_index, monoms, coeffs, params)
    res = _isgroebner1(ring, monoms, coeffs, params)
    res
end

# isgroebner for Finite fields
@timeit function _isgroebner1(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params
) where {M <: Monom, C <: CoeffZp}
    basis, pairset, hashtable = f4_initialize_structs(ring, monoms, coeffs, params)
    res = f4_isgroebner!(ring, basis, pairset, hashtable, params.arithmetic)
    res
end

# isgroebner for Rational numbers
@timeit function _isgroebner1(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params
) where {M <: Monom, C <: CoeffQQ}
    buffer = CoefficientBuffer()
    basis, pairset, hashtable = f4_initialize_structs(ring, monoms, coeffs, params)
    # If an honest computation over the rationals is needed
    if params.certify_check
        @log :misc """
        Keyword argument `certify=true` was provided. 
        Checking that the given input is a Groebner basis directly over the rationals"""
        flag = f4_isgroebner!(ring, basis, pairset, hashtable, params.arithmetic)
        return flag
    end
    # Otherwise, check modulo a prime
    @log :misc "Checking if a Grobner basis modulo a prime"
    buffer = CoefficientBuffer()
    @log :misc "Clearning denominators in the input generators"
    basis_zz = clear_denominators!(buffer, basis, deepcopy=false)
    luckyprimes = LuckyPrimes(basis_zz.coeffs)
    prime = next_check_prime!(luckyprimes)
    @log :misc "Reducing input generators modulo $prime"
    ring_ff, basis_ff = reduce_modulo_p!(buffer, ring, basis_zz, prime, deepcopy=true)
    arithmetic = select_arithmetic(CoeffModular, prime, :auto, false)
    flag = f4_isgroebner!(ring_ff, basis_ff, pairset, hashtable, arithmetic)
    flag
end
