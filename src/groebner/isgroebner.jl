
function isgroebner(
    polynomials;
    ordering=InputOrdering(),
    certify=false,
    seed=42,
    loglevel=Logging.Warn
)

    #= set the logger =#
    prev_logger = Logging.global_logger(ConsoleLogger(stderr, loglevel))

    #= extract ring information, exponents and coefficients
       from input polynomials =#
    # Copies input, so that polys would not be changed itself.
    ring, exps, coeffs =
        convert_to_internal(default_safe_representation(), polynomials, ordering)

    #= check and set algorithm parameters =#
    metainfo = set_metaparameters(ring, ordering, certify, :exact, rng)
    # now ring stores computation ordering
    # metainfo is now a struct to store target ordering

    iszerobasis = remove_zeros_from_input!(ring, exps, coeffs)
    iszerobasis && (return true)

    #= change input ordering if needed =#
    newring = assure_ordering!(ring, exps, coeffs, metainfo.targetord)

    #= check if groebner basis =#
    flag = isgroebner(newring, exps, coeffs, metainfo)

    #=
    Assuming ordering of `bexps` here matches `ring.ord`
    =#

    #= revert logger =#
    Logging.global_logger(prev_logger)

    flag
end

function isgroebner_f4!(
    ring::PolyRing,
    basis::Basis{C},
    ht::MonomialHashtable{M}
) where {M, C <: Coeff}
    matrix = initialize_matrix(ring, C)
    symbol_ht = hashtable(ht)
    update_ht = hashtable(ht)

    pairset = initialize_pairset(powertype(M))

    plcm = Vector{MonomIdx}(undef, 0)
    update!(pairset, basis, ht, update_ht, plcm)

    if pairset.load == 0
        return true
    end

    select_normal!(pairset, basis, matrix, ht, symbol_ht, selectall=true)

    symbolic_preprocessing!(basis, matrix, ht, symbol_ht)

    column_to_monom_mapping!(matrix, symbol_ht)

    sort_matrix_upper_rows_decreasing!(matrix) # for pivots,  AB part
    sort_matrix_lower_rows_increasing!(matrix)

    exact_sparse_rref_isgroebner!(ring, matrix, basis)
end

function isgroebner_f4(
    ring::PolyRing,
    exps::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    rng
) where {M, C <: Coeff}
    basis, ht = initialize_structs(ring, exps, coeffs, rng, 2^16)

    isgroebner_f4!(ring, basis, ht)
end

#------------------------------------------------------------------------------
# Finite field isgroebner

function isgroebner(
    ring::PolyRing,
    exps::Vector{Vector{M}},
    coeffs::Vector{Vector{T}},
    meta
) where {M, T <: CoeffFF}
    isgroebner_f4(ring, exps, coeffs, meta.rng)
end

#------------------------------------------------------------------------------
# Rational numbers groebner

function isgroebner(
    ring::PolyRing,
    exps::Vector{Vector{M}},
    coeffs::Vector{Vector{T}},
    meta
) where {M, T <: CoeffQQ}
    # if randomized result is ok
    if !meta.guaranteedcheck
        coeffbuffer = CoeffBuffer()
        coeffs_zz = scale_denominators(coeffbuffer, coeffs)
        primes = PrimeTracker(coeffs_zz)
        goodprime = nextgoodprime!(primes)
        coeffs_ff = reduce_modulo(coeffbuffer, coeffs_zz, goodprime)
        ring.ch = goodprime
        isgroebner_f4(ring, exps, coeffs_ff, meta.rng)
    else # if proved result is needed, then compute in rationals
        isgroebner_f4(ring, exps, coeffs, meta.rng)
    end
end
