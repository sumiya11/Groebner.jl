
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
    basis, ht = initialize_structures(ring, basisexps, basiscoeffs, rng, tablesize)

    tobereduced, ht =
        initialize_structures_nf(ring, tobereducedexps, tobereducedcfs, rng, tablesize, ht)

    tobereduced = normal_form_f4!(ring, basis, ht, tobereduced)

    export_basis_data(tobereduced, ht)
end
