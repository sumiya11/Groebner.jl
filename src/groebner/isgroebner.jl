
function _isgroebner(polynomials, kws::KeywordsHandler)
    polynomial_repr =
        select_polynomial_representation(polynomials, kws, hint=:large_exponents)
    ring, monoms, coeffs = convert_to_internal(polynomial_repr, polynomials, kws)
    if isempty(monoms)
        @log "Input consisting of zero polynomials, which is a Groebner basis by our convention"
        return true
    end
    params = AlgorithmParameters(ring, kws)
    ring = change_ordering_if_needed!(ring, monoms, coeffs, params)
    _isgroebner(ring, monoms, coeffs, params)
end

# isgroebner for Finite fields
function _isgroebner(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params
) where {M, C <: CoeffFF}
    basis, pairset, hashtable = initialize_structs(ring, monoms, coeffs, params)
    f4_isgroebner!(ring, basis, pairset, hashtable)
end

# isgroebner for Rational numbers
function _isgroebner(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params
) where {M, C <: CoeffQQ}
    buffer = CoefficientBuffer()
    basis, pairset, hashtable = initialize_structs(ring, monoms, coeffs, params)
    # If an honest computation over the rationals is needed
    if params.guaranteed_check
        @log "Checking if a Grobner basis directly over the rationals"
        return f4_isgroebner!(ring, basis, pairset, hashtable)
    end
    # Otherwise, check modulo a prime
    @log "Checking if a Grobner basis modulo a prime"
    buffer = CoefficientBuffer()
    @log "Clearning denominators in the input generators"
    basis_zz = clear_denominators!(buffer, basis, deepcopy=false)
    luckyprimes = LuckyPrimes(basis_zz.coeffs)
    prime = nextgoodprime!(luckyprimes)
    @log "Reducing input generators modulo $prime"
    ring_ff, basis_ff = reduce_modulo_p!(buffer, ring, basis_zz, prime, deepcopy=true)
    f4_isgroebner!(ring_ff, basis_ff, pairset, hashtable)
end
