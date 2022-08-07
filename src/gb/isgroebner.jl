
#------------------------------------------------------------------------------

function clean_input_isgroebner!(ring::PolyRing,
    exps::Vector{Vector{ExponentVector}},
    coeffs::Vector{CoeffsVector{T}}) where {T}
    clean_input_groebner!(ring, exps, coeffs)
end

#------------------------------------------------------------------------------

function isgroebner_f4!(ring::PolyRing,
    basis::Basis{C},
    ht::MonomialHashtable) where {C<:Coeff}

    matrix = initialize_matrix(ring, C)
    symbol_ht = initialize_secondary_hash_table(ht)
    update_ht = initialize_secondary_hash_table(ht)

    pairset = initialize_pairset()

    plcm = Vector{ExponentIdx}(undef, 0)
    update!(pairset, basis, ht, update_ht, plcm)

    if pairset.load == 0
        return true
    end

    select_isgroebner!(pairset, basis, matrix, ht, symbol_ht)

    symbolic_preprocessing!(basis, matrix, ht, symbol_ht)

    convert_hashes_to_columns!(matrix, symbol_ht)

    sort_matrix_rows_increasing!(matrix)
    sort_matrix_rows_decreasing!(matrix) # for pivots,  AB part

    exact_sparse_linear_algebra_isgroebner!(ring, matrix, basis)
end

function isgroebner_f4(
    ring::PolyRing,
    exps::Vector{Vector{ExponentVector}},
    coeffs::Vector{Vector{C}},
    rng) where {C<:Coeff}

    basis, ht = initialize_structures(ring, exps,
        coeffs, rng, 2^16)

    isgroebner_f4!(ring, basis, ht)
end

#######################################
# Finite field isgroebner

function isgroebner(
        ring::PolyRing,
        exps::Vector{Vector{ExponentVector}},
        coeffs::Vector{Vector{T}},
        meta) where {T<:CoeffFF}

    isgroebner_f4(ring, exps, coeffs, meta.rng)
end

#######################################
# Rational field groebner

function isgroebner(
        ring::PolyRing,
        exps::Vector{Vector{ExponentVector}},
        coeffs::Vector{Vector{T}},
        meta) where {T<:CoeffQQ}
        # if randomized result is ok
    if !meta.guaranteedcheck
        coeffbuffer  = CoeffBuffer()
        coeffs_zz = scale_denominators(coeffbuffer, coeffs)
        primes = PrimeTracker(coeffs_zz)
        goodprime = nextgoodprime!(primes)
        coeffs_ff = reduce_modulo(coeffbuffer, coeffs_zz, goodprime)
        ring.ch = goodprime
        isgroebner_f4(ring, exps, coeffs_ff, meta.rng)
    else # if proved result is needed, compute in rationals
        isgroebner_f4(ring, exps, coeffs, meta.rng)
    end
end