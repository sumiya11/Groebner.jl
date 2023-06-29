
function _normalform(basis, tobereduced, kws::KeywordsHandler) end

function normalform(basis::AbstractVector, tobereduced; kws...)
    iszero(tobereduced) && return tobereduced
    first(normalform(basis, [tobereduced], kws))
end

function normalform(basis::AbstractVector, tobereduced::AbstractVector, kws::KeywordsHandler)
    #= set the logger =#
    prev_logger = Logging.global_logger(ConsoleLogger(stderr, loglevel))

    check && _check_isgroebner(basispolys)

    #= extract ring information, exponents and coefficients
       from input basis polynomials =#
    # Copies input, so that polys would not be changed itself.
    ring1, basisexps, basiscoeffs =
        convert_to_internal(default_safe_representation(), basispolys, ordering)
    ring2, tbrexps, tbrcoeffs =
        convert_to_internal(default_safe_representation(), tobereduced, ordering)

    @assert ring1.nvars == ring2.nvars && ring1.ch == ring2.ch
    @assert ring1.ord == ring2.ord

    ring = ring1

    #= check and set algorithm parameters =#
    metainfo = set_metaparameters(ring, ordering, false, :exact, rng)

    iszerobasis = remove_zeros_from_input!(ring, basisexps, basiscoeffs)
    iszerobasis &&
        (return convert_to_output(ring, tobereduced, tbrexps, tbrcoeffs, metainfo))

    #= change input ordering if needed =#
    newring = assure_ordering!(ring, basisexps, basiscoeffs, metainfo.targetord)
    newring = assure_ordering!(ring, tbrexps, tbrcoeffs, metainfo.targetord)

    # We assume basispolys is already a Groebner basis! #

    #= compute the groebner basis =#
    bexps, bcoeffs =
        normal_form_f4(newring, basisexps, basiscoeffs, tbrexps, tbrcoeffs, rng)

    #=
    Assuming ordering of `bexps` here matches `newring.ord`
    =#

    #= revert logger =#
    Logging.global_logger(prev_logger)

    # ring contains ordering of computation, it is the requested ordering
    #= convert result back to representation of input =#
    convert_to_output(newring, tobereduced, bexps, bcoeffs, metainfo)
end

function normal_form_f4!(
    ring::PolyRing,
    basis::Basis{C},
    ht::MonomialHashtable,
    tobereduced::Basis{C}
) where {C <: Coeff}
    matrix = initialize_matrix(ring, C)
    symbol_ht = initialize_secondary_hash_table(ht)

    select_tobereduced!(basis, tobereduced, matrix, symbol_ht, ht)

    symbolic_preprocessing!(basis, matrix, ht, symbol_ht)

    column_to_monom_mapping!(matrix, symbol_ht)

    sort_matrix_upper_rows_decreasing!(matrix) # for pivots,  AB part

    exact_sparse_rref_nf!(ring, matrix, tobereduced, basis)

    convert_rows_to_basis_elements_nf!(matrix, tobereduced, ht, symbol_ht)

    tobereduced
end

function normal_form_f4(
    ring::PolyRing,
    basisexps::Vector{Vector{M}},
    basiscoeffs::Vector{Vector{C}},
    tobereducedexps::Vector{Vector{M}},
    tobereducedcfs::Vector{Vector{C}},
    rng::Rng,
    tablesize::Int=2^16
) where {M, Rng, C <: Coeff}
    basis, ht = initialize_structs(ring, basisexps, basiscoeffs, rng, tablesize)

    tobereduced, ht =
        initialize_structs_nf(ring, tobereducedexps, tobereducedcfs, rng, tablesize, ht)

    tobereduced = normal_form_f4!(ring, basis, ht, tobereduced)

    export_basis_data(tobereduced, ht)
end
